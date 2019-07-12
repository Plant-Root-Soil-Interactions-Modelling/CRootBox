import unittest
import py_rootbox as rb


class TestOrgan(unittest.TestCase):

    def hand_example(self):
        """ an example used in the tests below, a hand with two fingers """
        self.human1 = rb.Organism()  # same example as in test_constructor ...
        otp = rb.OrganTypeParameter(self.human1)
        self.human1.setOrganTypeParameter(otp)
        op = otp.realize()
        self.hand = rb.Organ(self.human1.getOrganIndex(), op, True, True, 0, 15., False, 0)
        self.hand.setOrganism(self.human1)
        self.thumb = rb.Organ(self.human1, self.hand, 0, 0, 4)  # delayedfor 4 days
        self.little_finger = rb.Organ(self.human1, self.hand, 0, 0, 3)  # delayed for 3 days
        self.hand.addChild(self.thumb)
        self.hand.addChild(self.little_finger)

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

    def test_constructors(self):
        """ tests three different kinds of constructors """
        human1 = rb.Organism()
        otp = rb.OrganTypeParameter(human1)
        human1.setOrganTypeParameter(otp)
        op = otp.realize()
        # 1. constructor from scratch
        hand = rb.Organ(human1.getOrganIndex(), op, True, True, 0, 15., False, 0)
        hand.setOrganism(human1)
        # 2. used in simulation (must have parent, since there is no nullptr in Pyhton)
        thumb = rb.Organ(human1, hand, 0, 0, 4)
        little_finger = rb.Organ(human1, hand, 0, 0, 3)
        hand.addChild(thumb)
        hand.addChild(little_finger)
        # 3. deep copy (with a factory function)
        human2 = rb.Organism()
        hand2 = hand.copy(human2)
        self.assertEqual(str(hand), str(hand2), "deep copy: the organs shold be equal")
        self.assertIsNot(hand.getParam(), hand2.getParam(), "deep copy: organs have same parameter set")

    def test_simulation(self):
        """ tests if the ages agree after a simulate call of 10 days"""
        self.hand_example()
        self.hand.simulate(10)
        self.assertEqual(self.hand.getAge(), 10, "wrong organ age")
        self.assertEqual(self.thumb.getAge(), 6, "wrong organ age")
        self.assertEqual(self.little_finger.getAge(), 7, "wrong organ age")

    def test_sequential(self):
        """ tests if the the organ tree can be represented in a seqential list"""
        self.hand_example()
        self.add_nodes()  # only organs with number of nodes > 1 are considered
        ring = rb.Organ(self.human1, self.thumb, 0, 0, 4)  # add a ring to the thumb
        self.thumb.addChild(ring)
        ring.addNode(rb.Vector3d(0, -1, 1.6), self.thumb.getNodeId(1), 4)
        ring.addNode(rb.Vector3d(0, -1, 1.6), 4)
        organs = self.hand.getOrgans()
        self.assertEqual(len(organs), 4, "wrong number of organs")

    def test_geometry(self):
        """ tests if nodes can be retrieved from the organ """
        self.hand_example()
        self.add_nodes()
        organs = [self.hand, self.thumb, self.little_finger]
        non = [4, 2, 2]  # number of nodes in the hand example
        for oi, o in enumerate(organs):
            # print("Organ #", o.getId())
            self.assertEqual(o.getNumberOfNodes(), non[oi], "wrong number of nodes")
            self.assertEqual(len(o.getSegments()), non[oi] - 1, "wrong number of segments")
            for i in range(0, o.getNumberOfNodes()):
                # print("node", o.getNode(i), "id", o.getNodeId(i), "creation time", o.getNodeCT(i))
                self.assertEqual(o.getNode(i).x, 0, "x coordinate should be zero")

    def test_parameter(self):
        """ tests if parameters per organ can be accessed """
        self.hand_example()
        self.add_nodes()
        a1 = self.little_finger.getParameter("age")
        a2 = self.thumb.getParameter("age")
        self.assertEqual(a1, -3, "wrong age")
        self.assertEqual(a2, -4, "wrong age")
        o0 = self.hand.getParameter("order")
        o1 = self.little_finger.getParameter("order")
        o2 = self.thumb.getParameter("order")
        self.assertEqual(o0, 0, "wrong order")
        self.assertEqual(o1, 1, "wrong order")
        self.assertEqual(o2, 1, "wrong order")

    def test_dynamics(self):
        """ tests if nodes created in last time step are correct """  #
        self.hand_example()
        self.hand.simulate(1)
        self.add_nodes()
        n0 = self.hand.getNumberOfNodes() - self.hand.getOldNumberOfNodes()
        n1 = self.little_finger.getNumberOfNodes() - self.little_finger.getOldNumberOfNodes()
        n2 = self.thumb.getNumberOfNodes() - self.thumb.getOldNumberOfNodes()
        self.assertEqual(n0, 4, "wrong number of new nodes")
        self.assertEqual(n1, 2, "wrong number of new nodes")
        self.assertEqual(n2, 2, "wrong number of new nodes")
        self.hand.simulate(1)
        n0 = self.hand.getNumberOfNodes() - self.hand.getOldNumberOfNodes()
        n1 = self.little_finger.getNumberOfNodes() - self.little_finger.getOldNumberOfNodes()
        n2 = self.thumb.getNumberOfNodes() - self.thumb.getOldNumberOfNodes()
        self.assertEqual(n0, 0, "wrong number of new nodes")
        self.assertEqual(n1, 0, "wrong number of new nodes")
        self.assertEqual(n2, 0, "wrong number of new nodes")
        self.hand.simulate(1)
        self.little_finger.addNode(rb.Vector3d(0, 1, 1.6), 6)
        n0 = self.hand.getNumberOfNodes() - self.hand.getOldNumberOfNodes()
        n1 = self.little_finger.getNumberOfNodes() - self.little_finger.getOldNumberOfNodes()
        n2 = self.thumb.getNumberOfNodes() - self.thumb.getOldNumberOfNodes()
        self.assertEqual(n0, 0, "wrong number of new nodes")
        self.assertEqual(n1, 1, "wrong number of new nodes")
        self.assertEqual(n2, 0, "wrong number of new nodes")


if __name__ == '__main__':
    unittest.main()
