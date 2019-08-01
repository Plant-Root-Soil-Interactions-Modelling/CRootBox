import unittest
import py_rootbox as rb
from rsml import *
from rb_tools import *


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


class TestRootSystem(unittest.TestCase):

    def rs_example_rtp(self):
        """ an example used in some of the tests below, 100 basals with laterals """
        self.rs = rb.RootSystem()
        maxB, firstB, delayB = 100, 10., 3
        rsp = rb.RootSystemParameter()
        rsp.set(3., firstB, delayB, maxB, 0, 1.e9, 1.e9, 1.e9, 0., 0.)
        self.rs.setRootSystemParameter(rsp)
        p0 = rb.RootTypeParameter(self.rs)
        p0.name, p0.type, p0.la, p0.nob, p0.ln, p0.r, p0.dx = "taproot", 1, 10, 20, 89. / 19., 1, 0.5
        p0.lb = 2
        p0.successor = a2i([2])  # to rb.std_int_double_()
        p0.successorP = a2v([1.])  # rb.std_vector_double_()
        p1 = rb.RootTypeParameter(self.rs)
        p1.name, p1.type, p1.la, p1.ln, p1.r, p1.dx = "lateral", 2, 25, 0, 2, 0.1
        self.p0, self.p1, self.rsp = p0, p1, rsp  # Python will garbage collect them away, if not stored
        self.rs.setOrganTypeParameter(self.p0)  # the organism manages the type parameters
        self.rs.setOrganTypeParameter(self.p1)

    def rs_length_test(self, dt, l, subDt):
        """ simulates a root system and checks basal lengths against its analytic lengths @param l at times @param t"""
        self.rs.initialize()
        nl = []
        for t in dt:
            for i in range(0, subDt):
                self.rs.simulate(t / subDt)
            ll = v2a(self.rs.getParameter("length"))
            types = v2a(self.rs.getParameter("type"))
            sl = 0  # summed length of basal roots
            for i, l_ in enumerate(ll):
                if (types[i] == 4):  # basal type
                    sl += l_
            nl.append(float(sl))
        for i in range(0, len(dt)):  # Check lengthes
            self.assertAlmostEqual(l[i], nl[i], 10, "numeric and analytic lengths do not agree (parameter length)")

    def rs_ct_test(self, dt, ct, subDt):
        """ simulates a root system and checks creation times analytic @param ct at times @param t"""
        self.rs.initialize()
        nl = []
        for t in dt:
            for i in range(0, subDt):
                self.rs.simulate(t / subDt)
            # creation times
            cts1 = v2a(self.rs.getParameter("creationTime"))
            poly_ct = self.rs.getPolylineCTs()
            cts2 = []
            for p in poly_ct:
                cts2.append(p[0])
            self.assertEqual(cts1.shape[0], len(cts2), "creation times: sizes are wrong")
            for i in range(0, len(cts2)):
                self.assertAlmostEqual(float(cts1[i]), float(cts2[i]), 10,
                                       "creation times: numeric creation times of polylines and parameter do not agree")
            types = v2a(self.rs.getParameter("type"))
            basal_ct = []
            for i, ct_ in enumerate(cts1):
                if types[i] == 4:
                    basal_ct.append(float(ct_))
            basal_ct.sort()
            for i in range(0, len(basal_ct)):
                self.assertAlmostEqual(float(basal_ct[i]), float(ct[i]), 10,
                                       "creation times: numeric and analytic creation times for basal roots")
            # tip times. as long the tips are active, tip creation time equals simulation time
            poly_ct = self.rs.getPolylineCTs()
            cts2 = []
            for p in poly_ct:
                cts2.append(p[-1])
            simtime = self.rs.getSimTime()
            for i in range(0, len(cts2)):
                self.assertAlmostEqual(float(cts2[i]), simtime, 10, "creation times: tip has wrong creation time")

    def test_length_no_laterals(self):
        """ run a simulation with a fibrous root system and compares to analytic lengths"""
        self.rs_example_rtp()
        times = np.array([0., 7., 15., 30., 60.])
        dt = np.diff(times)
        times = times[1:]
        etB = np.array(range(self.rsp.maxB)) * self.rsp.delayB + np.ones(self.rsp.maxB) * self.rsp.firstB  # basal root emergence times
        bl = np.zeros(times.size)  # summed root lengths
        for j, t in enumerate(times):
            i = 0  # basal root counter
            while t - etB[i] > 0:
                bl[j] += rootLength(t - etB[i], self.p0.r, self.p0.getK())
                i += 1
        self.rs_length_test(dt, bl, 1)
        self.rs_length_test(dt, bl, 100)

    def test_times_no_laterals(self):
        """ run a simulation with a fibrous root system and checks creation times """
        self.rs_example_rtp()
        times = np.array([0., 7., 15., 30., 60.])
        dt = np.diff(times)
        times = times[1:]
        ctB = np.array(range(self.rsp.maxB)) * self.rsp.delayB + np.ones(self.rsp.maxB) * self.rsp.firstB  # basal root emergence times
        # numeric solution
        self.rs_ct_test(dt, ctB, 1)
        self.rs_ct_test(dt, ctB, 100)

    def test_length_with_laterals(self):
        """ run a simulation with a fibrous root system and compares to analytic lengths"""
#         name = "Anagallis_femina_Leitner_2010"
#         rs = rb.RootSystem()
#         rs.openFile(name)
#         rs.initialize()
#         rs.simulate(7)
#         # todo
        pass
    # TODO

    def test_copy(self):
        """ checks if the root system can be copied, and if randomness works """
        seed = 100
        name = "Brassica_oleracea_Vansteenkiste_2014"
        rs = rb.RootSystem()  # the original
        rs.openFile(name)
        rs.setSeed(seed)
        rs.initialize()
        rs2 = rb.RootSystem(rs)  # copy root system
        self.assertIsNot(rs2, rs, "copy: not a copy")
        n1 = rs.rand()
        self.assertEqual(rs2.rand(), n1, "copy: random generator seed was not copied")
        rs.simulate(10)
        rs2.simulate(10)
        n2 = rs.rand()
        self.assertEqual(rs2.rand(), n2, "copy: simulation not deterministic")
        rs3 = rb.RootSystem()  # rebuild same
        rs3.openFile(name)
        rs3.setSeed(seed)
        rs3.initialize()
        self.assertEqual(rs3.rand(), n1, "copy: random generator seed was not copied")
        rs3.simulate(10)
        self.assertEqual(rs3.rand(), n2, "copy: simulation not deterministic")

    def test_polylines(self):
        """checks if the polylines have the right tips and bases """
        name = "Brassica_napus_a_Leitner_2010"
        rs = rb.RootSystem()
        rs.openFile(name)
        rs.initialize()
        rs.simulate(7)  # days young
        polylines = rs.getPolylines()  # Use polyline representation of the roots
        bases = np.zeros((len(polylines), 3))
        tips = np.zeros((len(polylines), 3))
        for i, r in enumerate(polylines):
            bases[i, :] = [r[0].x, r[0].y, r[0].z]
            tips[i, :] = [r[-1].x, r[-1].y, r[-1].z]
        nodes = vv2a(rs.getNodes())  # Or, use node indices to find tip or base nodes
        tipI = rs.getRootTips()
        baseI = rs.getRootBases()
        uneq = np.sum(nodes[baseI, :] != bases) + np.sum(nodes[tipI, :] != tips)
        self.assertEqual(uneq, 0, "polylines: tips or base nodes do not agree")

    def test_root_type_parameters(self):
        """ root type parameters xml read and write """
        self.rs_example_rtp()
#         print(self.p0.__str__(False))
#         print(self.p1.__str__(False))
        print(rb.Organism.organTypeName(self.p0.organType))
        self.rs.writeParameters("test_parameters.xml", "RootBox", False)  # include comments
        rs1 = rb.RootSystem
        # rs1.readParameters("test_parameters.xml", "RootBox")
        # TODO

#
#     def test_adjacency_matrix(self):
#         """ builds an adjacency matrix, and checks if everything is connected"""
#         pass
#
    def test_dynamics(self):
        """ incremental root system growth like needed for coupling"""
        name = "Anagallis_femina_Leitner_2010"  # "maize_p2"  # "Anagallis_femina_Leitner_2010"  # "Zea_mays_4_Leitner_2014"
        rs = rb.RootSystem()
        rs.openFile(name)
        rs.initialize()
        simtime = 60  # days
        dt = 1
        N = round(simtime / dt)
        nodes = vv2a(rs.getNodes())  # contains the initial nodes of tap, basal and shootborne roots
        seg = np.array([], dtype = np.int64).reshape(0, 2)
        cts = v2a(rs.getSegmentCTs())
        nonm = 0
        for i in range(0, N):
            rs.simulate(dt, False)
            uni = v2ai(rs.getUpdatedNodeIndices())  # MOVED NODES
            unodes = vv2a(rs.getUpdatedNodes())
            nodes[uni] = unodes  # do the update
            nonm += uni.shape[0]
            newnodes2 = rs.getNewNodes()  # NEW NODES
            newnodes = vv2a(newnodes2)
            newsegs = seg2a(rs.getNewSegments())  # NEW SEGS
            newcts = v2a(rs.getNewSegmentCTs())
            if len(newnodes) != 0:
                nodes = np.vstack((nodes, newnodes))
            if len(newsegs) != 0:
                seg = np.vstack((seg, newsegs))
                cts = np.vstack((cts, newcts))

        nodes_ = vv2a(rs.getNodes())
        nodeCTs_ = v2a(rs.getNodeCTs())
        seg_ = seg2a(rs.getSegments())
        sct_ = v2a(rs.getSegmentCTs())
#         print("Creation times range from ", np.min(cts), " to ", np.max(cts), "len", cts.shape)
#         print("Creation times range from ", np.min(sct_), " to ", np.max(sct_), "len", sct_.shape)
#         print("seg:", seg_[0])
#         print(nodes_[seg[0][0]], nodes_[seg[0][1]])
#         print(nodeCTs_[seg[0][0]], nodeCTs_[seg[0][1]])
#         print("---")
        self.assertEqual(nodes_.shape, nodes.shape, "incremental growth: node lists are not equal")
        ind = np.argwhere(nodes_[:, 1] != nodes[:, 1])
        for i in ind:
            print(i, nodes_[i], "!=", nodes[i])
        uneq = np.sum(nodes_ != nodes) / 3
        self.assertEqual(uneq, 0, "incremental growth: node lists are not equal")
        self.assertEqual(seg_.shape, seg.shape, "incremental growth: segment lists are not equal")
        seg = np.sort(seg, axis = 0)  # per default along the last axis
        seg_ = np.sort(seg_, axis = 0)
        uneq = np.sum(seg_ != seg) / 2
        self.assertEqual(uneq, 0, "incremental growth: segment lists are not equal")

    def test_rsml(self):
        """ checks rsml functionality with Python rsml reader """
        name = "Anagallis_femina_Leitner_2010"
        rs = rb.RootSystem()
        rs.openFile(name)
        rs.initialize()
        simtime = 60
        rs.simulate(simtime)
        rs.writeRSML(name + ".rsml")
        pl, props, funcs = read_rsml(name + ".rsml")
        # todo

#     def test_stack(self):
#         """ checks if push and pop are working """


if __name__ == '__main__':
    unittest.main()
