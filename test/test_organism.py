import unittest
import py_rootbox as rb
from rsml import *


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
        self.human1.getNodeIndex()
        self.human1.getNodeIndex()  # to make global and organ index disagree
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
        self.hand_example()
        human2 = rb.Organism(self.human1)  # copy constructor
        self.assertIsNot(self.human1, human2, "copy: not a copy")
        self.assertEqual(self.human1.rand(), human2.rand(), "copy: random generator seed was not copied")
        # todo check organs
        # check otps

    def test_organ_type_parameters(self):
        pass

    def test_simulation(self):
        pass

    def test_geometry(self):
        pass

    def test_dynamics(self):  #
        pass

    def test_rsml(self):
        """ checks rmsl functionality with Python rsml reader
        """
        self.hand_example()
        self.add_nodes()
        self.human1.addOrgan(self.hand)
        self.human1.writeRSML("organism.rsml")
        pl, props, funcs = read_rsml("organism.rsml")
        pl2 = [[[0.0, 0.0, 0.0], [0.0, 0.0, 1.5], [0.0, -1.0, 1.6], [0.0, 1.0, 1.6]], [[0.0, -1.0, 1.6], [0.0, -2.0, 2.5]],
               [[0.0, 1.0, 1.6], [0.0, 1.7, 2.5]], [[0.0, 0.0, 0.0], [0.0, 0.0, 1.5], [0.0, -1.0, 1.6], [0.0, 1.0, 1.6]],
               [[0.0, -1.0, 1.6], [0.0, -2.0, 2.5]], [[0.0, 1.0, 1.6], [0.0, 1.7, 2.5]]]
        self.assertEqual(pl, pl2, "rsml: polylines are not equal")
        # todo test everything


if __name__ == '__main__':
    unittest.main()
