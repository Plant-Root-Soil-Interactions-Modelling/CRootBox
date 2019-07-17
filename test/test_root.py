import unittest
import py_rootbox as rb
import numpy as np


def rootAge(l, r, k):  # root age at a certain length
    return -np.log(1 - l / k) * k / r


def rootLength(t, r, k):  # root length at a certain age
    return k * (1 - np.exp(-r * t / k))


def rootLateralLength(t, et, r, k):  # length of first order laterals (without second order laterals)
    i, l = 0, 0
    while et[i] < t:
        age = t - et[i]
        l += rootLength(age, r, k)
        i += 1
    return l


class TestRoot(unittest.TestCase):

    def root_example_rtp(self):
        """ an example used in the tests below, a main root with laterals """
        self.plant = rb.Organism()  # Root has no dependency on RootSystem anymore
        p0 = rb.RootTypeParameter(self.plant)
        p1 = rb.RootTypeParameter(self.plant)
        p0.name = "taproot"
        p0.type = 1
        p0.lb = 1
        p0.la = 10
        p0.nob = 20
        p0.ln = 89. / 19.
        p0.r = 1
        p0.dx = 0.5
        l = rb.std_vector_int_()  # set up successors
        l.append(2)
        p0.successor = l
        l = rb.std_vector_double_()
        l.append(1)
        p0.successorP = l
        p1.name = "lateral"
        p1.type = 2
        p1.la = 25
        p1.ln = 0
        p1.r = 2
        p1.dx = 0.1
        self.p0 = p0
        self.p1 = p1
        self.plant.setOrganTypeParameter(self.p0)  # the organism manages the type parameters
        self.plant.setOrganTypeParameter(self.p1)
        self.param0 = self.p0.realize()  # set up root by hand (without a root system)
        self.param0 .la = 0  # its important parent has zero length, otherwise creation times are messed up
        self.param0 .lb = 0
        # param0 is stored, because otherwise garbage collection deletes it, an program will crash <---
        parentroot = rb.Root(1, self.param0, True, True, 0., 0., rb.Vector3d(0, 0, -1), 0, 0, False, 0)
        parentroot.setOrganism(self.plant)
        parentroot.addNode(rb.Vector3d(0, 0, -3), 0)  # there is no nullptr in Python
        self.root = rb.Root(self.plant, self.p0.subType, rb.Vector3d(0, 0, -1), 0, parentroot, 0, 0)
        self.root.setOrganism(self.plant)

    def root_length_test(self, dt, l, subDt):
        """ simulates a single root and checks length against analytic length """
        nl, nl2, non, meanDX = [], [], [], []
        for t in dt:
            for i in range(0, subDt):
                self.root.simulate(t / subDt)
            nl.append(self.root.getParameter("length"))
            non.append(self.root.getNumberOfNodes())
            meanDX.append(nl[-1] / non[-1])
            # length from geometry
            poly = np.zeros((non[-1], 3))  #
            for i in range(0, non[-1]):
                v = self.root.getNode(i)
                poly[i, 0] = v.x
                poly[i, 1] = v.y
                poly[i, 2] = v.z
            d = np.diff(poly, axis = 0)
            sd = np.sqrt((d ** 2).sum(axis = 1))
            nl2.append(sum(sd))
        for i in range(0, len(dt)):
            self.assertAlmostEqual(l[i], nl[i], 10, "numeric and analytic lengths do not agree in time step " + str(i + 1))
            self.assertAlmostEqual(l[i], nl2[i], 10, "numeric and analytic lengths do not agree in time step " + str(i + 1))
            self.assertLessEqual(meanDX[i], 0.5, "axial resolution dx is too large")
            self.assertLessEqual(0.25, meanDX[i], "axial resolution dx is unexpected small")

    def test_constructors(self):
        """ tests three different kinds of constructors """
        self.root_example_rtp()
        # 1. constructor from scratch
        param = self.p0.realize()
        root = rb.Root(1, param, True, True, 0., 0., rb.Vector3d(0, 0, -1), 0, 0, False, 0)
        root.setOrganism(self.plant)
        root.addNode(rb.Vector3d(0, 0, -3), 0)  # parent must have at least one nodes
        # 2. used in simulation (must have parent, since there is no nullptr in Pyhton)
        root2 = rb.Root(self.plant, self.p1.subType, rb.Vector3d(0, 0, -1), 0, root, 0, 0)
        root.addChild(root2);
        # 3. deep copy (with a factory function)
        plant2 = rb.Organism()
        root3 = root.copy(plant2)
        self.assertEqual(str(root), str(root3), "deep copy: the organs shold be equal")
        self.assertIsNot(root.getParam(), root3.getParam(), "deep copy: organs have same parameter set")
        # TODO check if OTP were copied

    def test_root_length(self):
        """ tests if numerical root length agrees with analytic solutions at 4 points in time with two scales of dt"""
        self.root_example_rtp()
        times = np.array([0., 7., 15., 30., 60.])
        dt = np.diff(times)
        k = self.root.param().getK()  # maximal root length
        self.assertAlmostEqual(k, 100, 12, "example root has wrong maximal length")
        l = rootLength(times[1:], self.p0.r, k)  # analytical root length
        root = self.root.copy(self.plant)
        self.root_length_test(dt, l, 1)  # large dt
        self.root = root
        self.root_length_test(dt, l, 1000)  # very fine dt

    def test_root_length_including_laterals(self):
        """ tests if numerical root length agrees with analytic solution including laterals """
        self.root_example_rtp()
        times = np.array([0., 7., 15., 30., 60.])
        dt = np.diff(times)
        p = self.root.param()  # rename
        k = p.getK()
        et = np.zeros((p.nob))
        l = 0
        et[0] = rootAge(p.la + p.lb + l, p.r, k)
        for i in range(0, p.nob - 1):  # calculate lateral emergence times
            l += p.ln[i]
            et[i + 1] = rootAge(p.la + p.lb + l, p.r, k + 1e-12)
        l = rootLength(times[1:], p.r, k)  # zero order lengths
        l1 = []
        r2 = self.p1.r
        k2 = self.p1.la  # consists of lateral zone only
        for t in times[1:]:
            l1.append(rootLateralLength(t, et, r2, k2))
        analytic_total = l + l1

        for subDX in [1, 1000]:
            numeric_total = []
            for t in times[1:]:
                root = self.root.copy(self.plant)
                self.root_length_test([t], [rootLength(t, p.r, k)], subDX)
                organs = self.root.getOrgans()
                nl = 0
                for o in organs:
                    nl += o.getParameter("length")
                numeric_total.append(nl);
                self.root = root
            for i in range(0, len(times[1:])):
                self.assertAlmostEqual(numeric_total[i], analytic_total[i], 10, "numeric and analytic total lengths do not agree in time step " + str(i + 1))

    def test_geometry(self):
        """ tests if nodes can be retrieved from the organ """
        # make plot for plausibility

    def test_parameter(self):
        """ tests some parameters on sequential organ list """
        self.root_example_rtp()
        self.root.simulate(30)
        organs = self.root.getOrgans()
        type, age, radius, order, ct = [], [], [], [], []
        for o in organs:
            type.append(o.getParameter("subType"))
            age.append(o.getParameter("age"))
            ct.append(o.getParameter("creationTime"))
            radius.append(o.getParameter("radius"))
            order.append(o.getParameter("order"))
        self.assertEqual(type, [1.0, 2.0, 2.0, 2.0, 2.0], "getParameter: unexpected root sub types")
        self.assertEqual(order, [1.0, 2.0, 2.0, 2.0, 2.0], "getParameter: unexpected root sub types")  # +1, because of artificial parent root
        for i in range(0, 5):
            self.assertEqual(age[i], 30 - ct[i], "getParameter: unexpected root sub types")  # +1, because of artificial parent root
#

    def test_dynamics(self):
        """ tests if nodes created in last time step are correct """  #
        self.root_example_rtp()
        r = self.root
        r.simulate(.5, True)
        self.assertEqual(r.hasMoved(), False, "dynamics: node movement during first step")
        r.simulate(1e-1, True)
        self.assertEqual(r.hasMoved(), False, "dynamics: movement, but previous node at axial resolution")
        r.simulate(1e-1, True)
        self.assertEqual(r.hasMoved(), True, "dynamics: node was expected to move, but did not")
        r.simulate(2.4, True)
        self.assertEqual(r.getNumberOfNodes() - r.getOldNumberOfNodes(), 6, "dynamics: unexcpected number of new nodes")


if __name__ == '__main__':
    unittest.main()
