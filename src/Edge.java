import java.util.ArrayList;

/**
 *  
 * Implementation of an Edge.
 * 
 * @see Graph.java
 * 
 */
public class Edge implements Comparable<Edge>{

	/**
	 * Edge attributes
	 */
	String node1;
	String node2;
	double weight; 
	double numActivities; 
	double originalNumActivities; //num activities may change during monte carlo simulations - keep original num activities for later
	ArrayList<Activity> activities;
	
	/**
	 * Constructs an edge
	 * @param node1
	 * @param node2
	 * @param weight
	 */
	public Edge(String node1, String node2, double weight) {
		this.node1 = node1;
		this.node2 = node2;
		this.weight = weight;
		numActivities = 0;
		this.originalNumActivities = numActivities;
		activities = new ArrayList<Activity>();
	}

	@Override
	public String toString() {
		return node1 + " -> " + node2 + "(" + weight + "," + numActivities + ")";
	}

	/**
	 * Compares this edge to another edge. Returns 1 if this edge's weight is greater
	 * that the given edges weight, -1 if lesser and 0 if the edge weights are equal. 
	 * @param edge - the edge that this edge is being compared with
	 * @return 1 if this edge's weight is greater that the given edges weight, 
	 * -1 if lesser and 0 if the edge weights are equal. 
	 */
	public int compareTo(Edge edge) {
		if (this.weight > edge.weight) {
			return 1;
		}
		else if (this.weight < edge.weight) {
			return -1;
		}
		else {
			return 0;
		}
	}
	
	@Override
	public boolean equals(Object object) {
		// two edge are the same or reversed
		if (object instanceof Edge) {
			Edge edge = (Edge) object;
			if (this.node1 == edge.node1 && this.node2 == edge.node2 && this.weight == edge.weight) {
				return true;
			}
			else if (this.node1 == edge.node2 && this.node2 == edge.node1 && this.weight == edge.weight) {
				return true;
			}
			else {
				return false;
			}			
		}
		return false;
	}

	@Override
	protected Edge clone() throws CloneNotSupportedException {
		Edge edge = new Edge(this.node1, this.node2, this.weight);
		edge.numActivities = this.numActivities;
		edge.activities.addAll(this.activities);
		return edge;
	}	
	
	/** 
	 * Returns the number of activities on dynamic path
	 * @return the number of activities on dynamic path
	 */
	public double calculateNumActivitiesDynamicForSPTree() {
		if (node2.startsWith("a.")) return 1;
		else return 0;
	}		
	
	
	
}
