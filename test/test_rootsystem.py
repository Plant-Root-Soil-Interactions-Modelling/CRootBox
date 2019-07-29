import unittest
import py_rootbox as rb
from rsml import *
from rb_tools import *


class TestRootSystem(unittest.TestCase):

    def root_example_rtp(self):
        """ an example used in the tests below, a main root with laterals """
        self.rs = rb.RootSystem()
        p0 = rb.RootTypeParameter(self.rs)
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
        p1 = rb.RootTypeParameter(self.rs)
        p1.name = "lateral"
        p1.type = 2
        p1.la = 25
        p1.ln = 0
        p1.r = 2
        p1.dx = 0.1
        self.p0, self.p1 = p0, p1  # Python will garbage collect them away, if not stored
        self.rs.setOrganTypeParameter(self.p0)  # the organism manages the type parameters
        self.rs.setOrganTypeParameter(self.p1)

    def root_length_test(self, dt, l, subDt):
        """ simulates a single root and checks length against its analytic length """
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
            self.assertAlmostEqual(l[i], nl[i], 10, "numeric and analytic lengths do not agree")
            self.assertAlmostEqual(l[i], nl2[i], 10, "numeric and analytic lengths do not agree")
            self.assertLessEqual(meanDX[i], 0.5, "axial resolution dx is too large")
            self.assertLessEqual(0.25, meanDX[i], "axial resolution dx is unexpected small")

    def test_root_type_parameters(self):
        """ root type parameters xml read and write """
        self.root_example_rtp()
        print(self.p0.__str__(False))
        print(self.p1.__str__(False))
        print(rb.Organism.organTypeName(self.p0.organType))
        self.rs.writeParameters("test_parameters.xml", "RootBox", False)  # include comments
        rs1 = rb.RootSystem
        # rs1.readParameters("test_parameters.xml", "RootBox")

#     def test_length_no_laterals(self):
#         """ run a simulation with a fibrous root system and compares to analytic lengths"""
#         name = "Anagallis_femina_Leitner_2010"
#         rs = rb.RootSystem()
#         rs.openFile(name)
#         rs.initialize()
#         rs.simulate(7)
#         # todo

#     def test_length_with_laterals(self):
#         """ run a simulation with a fibrous root system and compares to analytic lengths"""
#         name = "Anagallis_femina_Leitner_2010"
#         rs = rb.RootSystem()
#         rs.openFile(name)
#         rs.initialize()
#         rs.simulate(7)
#         # todo

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
        """checks if the polyliens have the right tips and bases """
        name = "Brassica_napus_a_Leitner_2010"
        rs = rb.RootSystem()
        rs.openFile(name)
        rs.initialize()
        rs.simulate(7)  # days young
        polylines = rs.getPolylines()  # Use polyline representation of the roots
        # todo getPolylinesCTs
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

#
#     def test_adjacency_matrix(self):
#         """ builds an adjacency matrix, and checks if everything is connected"""
#         pass
#
    def test_dynamics(self):
        """ incremental root system growth like needed for coupling"""

        # test currently does not run for Anagallis

        name = "Anagallis_femina_Leitner_2010"  # "Anagallis_femina_Leitner_2010"  # "Anagallis_femina_Leitner_2010"  # "Zea_mays_4_Leitner_2014"
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

#     def test_rsml(self):
#         """ checks rsml functionality with Python rsml reader """
#         name = "Anagallis_femina_Leitner_2010"
#         rs = rb.RootSystem()
#         rs.openFile(name)
#         rs.initialize()
#         simtime = 60
#         rs.simulate(simtime)
#         rs.writeRSML(name + ".rsml")
#         pl, props, funcs = read_rsml(name + ".rsml")
#         # todo

#     def test_stack(self):
#         """ checks if push and pop are working """


if __name__ == '__main__':
    unittest.main()
