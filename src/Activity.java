/**
 * An activity represents an object of interest in the network.
 */

public class Activity {

	/** 
	 * Attributes
	 */
	public String id;
	public Edge edge;
	public double x;
	public double y;
	public double distanceToEdgeNode1;
	public double distanceToEdgeNode2;
	
	public Path pathToNode1;
	public Path pathToNode2;
	
	public Activity (String id, double x, double y, Edge edge) {
		this.id = id;
		this.edge = edge;
		this.x = x;
		this.y = y;
		
		Node edgeNode1 = Graph.nodesInfoTable.get(edge.node1);
		double edgeNode1X = edgeNode1.x;
		double edgeNode1Y = edgeNode1.y;
		distanceToEdgeNode1 =  Math.sqrt(Math.pow((edgeNode1X-x),2)+Math.pow((edgeNode1Y-y),2));	
		
		Node edgeNode2 = Graph.nodesInfoTable.get(edge.node2);
		double edgeNode2X = edgeNode2.x;
		double edgeNode2Y = edgeNode2.y;
		distanceToEdgeNode2 =  Math.sqrt(Math.pow((edgeNode2X-x),2)+Math.pow((edgeNode2Y-y),2));	
	}		

	public void updateActivityInfo (double x, double y, Edge edge) {
		
		this.x = x;
		this.y = y;	
		this.edge = edge;
		
		Node edgeNode1 = Graph.nodesInfoTable.get(edge.node1);
		double edgeNode1X = edgeNode1.x;
		double edgeNode1Y = edgeNode1.y;
		distanceToEdgeNode1 =  Math.sqrt(Math.pow((edgeNode1X-x),2)+Math.pow((edgeNode1Y-y),2));	
		
		Node edgeNode2 = Graph.nodesInfoTable.get(edge.node2);
		double edgeNode2X = edgeNode2.x;
		double edgeNode2Y = edgeNode2.y;
		distanceToEdgeNode2 =  Math.sqrt(Math.pow((edgeNode2X-x),2)+Math.pow((edgeNode2Y-y),2));		
	}
	
	
	@Override	
	public String toString() {
		return "Activity " +id + " (" + x + ", " + y + ") is on edge " + edge.node1 + "->" + edge.node2 + " d=" + distanceToEdgeNode1 + "," + distanceToEdgeNode2 + "\n";
	}
	
}
