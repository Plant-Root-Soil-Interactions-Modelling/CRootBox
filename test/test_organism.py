import unittest
import py_rootbox as rb
import matplotlib.pyplot as plt
from rsml import *
from rb_tools import *


class TestOrganism(unittest.TestCase):

    def hand_example(self):
        """ an example used in the tests below (same as test_organ), a hand with two fingers """
        self.human1 = rb.Organism()  # same example as in test_constructor ...
        otp = rb.OrganTypeParameter(self.human1)
        self.human1.setOrganTypeParameter(otp)
        op = otp.realize()
        self.hand = rb.Organ(self.human1.getOrganIndex(), op, True, True, 0, 15., False, 0)
        self.hand.setOrganism(self.human1)
        self.thumb = rb.Organ(self.human1, self.hand, 0, 0, 4)  # delayed for 4 days
        self.little_finger = rb.Organ(self.human1, self.hand, 0, 0, 3)  # delayed for 3 days
        self.hand.addChild(self.thumb)
        self.hand.addChild(self.little_finger)
        self.human1.addOrgan(self.hand)

    def add_nodes(self):
        """ used in the tests below, adds nodes to the hand example """
        self.hand.addNode(rb.Vector3d(0, 0, 0), 0)
        self.hand.addNode(rb.Vector3d(0, 0, 1.5), 0)
        self.hand.addNode(rb.Vector3d(0, -1, 1.6), 0)  # thumb
        self.hand.addNode(rb.Vector3d(0, 1, 1.6), 0)  # little finger
        thumb = self.hand.getNodeId(2)
        lf = self.hand.getNodeId(3)
        self.thumb.addNode(rb.Vector3d(0, -1, 1.6), thumb, 4)
        self.thumb.addNode(rb.Vector3d(0, -2, 2.5), 4)
        self.little_finger.addNode(rb.Vector3d(0, 1, 1.6), lf, 3)
        self.little_finger.addNode(rb.Vector3d(0, 1.7, 2.5), 3)

    def test_copy(self):
        """ checks if the organism is properly copied """
        self.hand_example()
        self.add_nodes()
        human2 = rb.Organism(self.human1)  # copy constructor
        self.assertIsNot(self.human1, human2, "copy: not a copy")
        self.assertEqual(self.human1.rand(), human2.rand(), "copy: random generator seed was not copied")
        o1 = self.human1.getOrgans()  # check organs
        o2 = human2.getOrgans()
        self.assertEqual(len(o2), 3, "copy: unexpected number of organs")
        for i in range(0, len(o1)):
            self.assertIsNot(o1[i], o2[i], "copy: organ is not copied")
        p1 = self.human1.getOrganTypeParameter(0)
        p2 = human2.getOrganTypeParameter(0)
        for i in range(0, len(p1)):
            self.assertIsNot(p1[i], p2[i], "copy: OrganTypeParameters is not copied")

    def test_organ_type_parameters(self):
        """ test ability to set, get, read, and write type parameters """
        human1 = rb.Organism()  # same example as in test_constructor ...
        otp1 = rb.OrganTypeParameter(human1)
        otp1.name = "nose"
        otp1.subType = 1
        otp2 = rb.OrganTypeParameter(human1)
        otp2.subType = 2
        otp2.name = "eye"
        human1.setOrganTypeParameter(otp1)  # set
        human1.setOrganTypeParameter(otp2)
        otps = human1.getOrganTypeParameter(rb.OrganTypes.organ)
        self.assertEqual(otps[0].name, "nose", "otp: name not expected ")
        self.assertEqual(otps[1].name, "eye", "otp: name not expected ")
        otp3 = rb.OrganTypeParameter(human1)
        otp3.organType = rb.OrganTypes.root
        otp3.subType = 1
        otp3.name = "rootyhand"
        human1.setOrganTypeParameter(otp3)
        human1.writeParameters("human.xml")
        human2 = rb.Organism()  # read again
        prototype1 = rb.OrganTypeParameter(human2)
        prototype1.organType = rb.OrganTypes.organ
        prototype2 = rb.OrganTypeParameter(human2)
        prototype2.organType = rb.OrganTypes.root
        human2.setOrganTypeParameter(prototype1)  # set prototypes for reading, subTypes are overwritten if equal
        human2.setOrganTypeParameter(prototype2)
        human2.readParameters("human.xml")
        otp1 = human2.getOrganTypeParameter(rb.OrganTypes.organ, 1)
        otp2 = human2.getOrganTypeParameter(rb.OrganTypes.organ, 2)
        self.assertEqual(otp1.name, "nose", "otp: name not expected ")
        self.assertEqual(otp2.name, "eye", "otp: name not expected ")
        self.assertEqual(otp1.subType, 1, "otp: subType not expected ")
        self.assertEqual(otp2.subType, 2, "otp: subType not expected ")
        rtp = human2.getOrganTypeParameter(rb.OrganTypes.root, 1)
        self.assertEqual(rtp.name, "rootyhand", "otp: name not expected ")
        self.assertEqual(rtp.subType, 1, "otp: subType not expected ")

    def test_simulation(self):
        """ tests if the organs have the right age after simulation """
        self.hand_example()
        self.add_nodes()
        self.human1.simulate(365)  # happy birthday
        organs = self.human1.getOrgans()  # check organs
        self.assertEqual(organs[0].getAge(), 365, "simulation: age not expected ")
        self.assertEqual(organs[1].getAge(), 361, "simulation: age not expected ")
        self.assertEqual(organs[2].getAge(), 362, "simulation: age not expected ")
        self.human1.simulate(0.5)
        self.assertEqual(organs[0].getAge(), 365.5, "simulation: age not expected ")
        self.assertEqual(organs[1].getAge(), 361.5, "simulation: age not expected ")
        self.assertEqual(organs[2].getAge(), 362.5, "simulation: age not expected ")

    def test_geometry(self):
        """ tests ability to retrieve geometry """
        self.hand_example()
        self.add_nodes()
        nodes = vv2a(self.human1.getNodes())
        self.assertEqual(nodes.shape, (6, 3), "geometry: number of nodes unexpected")
        segs = seg2a(self.human1.getSegments())
        self.assertEqual(np.sum(np.sum(segs.flat != np.array([[0, 1], [1, 2], [2, 3], [2, 4], [3, 5]]).flat)), 0, "geometry: segments ids are unexcpected")

    def test_parameter(self):
        """ test if getParameter works """
        self.hand_example()
        self.add_nodes()
        organs = self.human1.getOrgans()
        age = v2a(self.human1.getParameter("age"))
        self.assertEqual(len(organs), age.shape[0] , "parameter: organ size unequal to parameter size")
        self.assertEqual(age[2], -3, "parameter: wrong ages")

    def test_rsml(self):
        """ checks rmsl functionality with Python rsml reader """
        self.hand_example()
        self.add_nodes()
        self.human1.writeRSML("organism.rsml")
        pl, props, funcs = read_rsml("organism.rsml")
        pl2 = [[[0.0, 0.0, 0.0], [0.0, 0.0, 1.5], [0.0, -1.0, 1.6], [0.0, 1.0, 1.6]], [[0.0, -2.0, 2.5]], [[0.0, 1.7, 2.5]]]
        self.assertEqual(pl, pl2, "rsml: polylines are not equal")
        self.assertEqual(props["age"], [0, -4, -3] , "rsml: polylines are not equal")


if __name__ == '__main__':
    unittest.main()
