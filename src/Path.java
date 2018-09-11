/**
 * Implementation of a path. <br><br>
 * 
 * A path is a sequence of connected edges. For example: 1->2->5 is a path consisting
 * of the edges 1->2 and 2->5.
 * 
 */

import java.math.BigDecimal;
import java.util.ArrayList;


public class Path implements Comparable<Path>{

	/**
	 * Attributes
	 */
	public ArrayList<Edge> edges;
	public int calculatedActivities;
	public double weight;
	public double numActivities;
	//public double numDynamicActivities;
	public double clusterActivities;	
	public double statistic;
	public double pVal;
	
	public double tempActivities;
	
	/**
	 * Constructs a path
	 */
	public Path() {
		edges = new ArrayList<Edge>();
		calculatedActivities = 0;
		weight = Double.MAX_VALUE;
		numActivities = 0;
		//numDynamicActivities = 0;
		clusterActivities = Double.MIN_VALUE;
		statistic = 0;
		pVal = 0;
	}
	
	/**
	 * Constructs a path
	 * @param edges - the edges that make up this path
	 */
	public Path(ArrayList<Edge> edges) {
		this.edges = edges;
		calculatedActivities = 0;
		weight = calculateWeight();
		numActivities = calculateNumActivities();
		//numDynamicActivities = 0;
		clusterActivities = Double.MIN_VALUE;
		statistic = 0;
		pVal = 0;
	}

	
	/**
	 * Returns the number of Activities on this path
	 * @return the number of Activities on this path
	 */
	public int getActivityCount(){
		int numActivities = 0;		
		for (int i = 0; i < edges.size(); i++) {
			numActivities += edges.get(i).numActivities;
		}	
		return numActivities;
	}		
	
	/** 
	 * Returns the weight of this path
	 * @return the weight of this path
	 */
	public double calculateWeight() {
		double weight = 0;
		for (int i = 0; i < edges.size(); i++){
			weight = weight + edges.get(i).weight;
		}
		return weight;
	}
	
	/** 
	 * Returns the number of activities on this path
	 * @return the number of activities on this path
	 */
	public double calculateNumActivities() {
		double numActivities = 0;
		for (int i = 0; i < edges.size(); i++){
			numActivities = numActivities + edges.get(i).numActivities;
		}
		return numActivities;
	}	
	
	/** 
	 * Returns the number of activities on dynamic path
	 * @return the number of activities on dynamic path
	 */
	public double calculateNumActivitiesDynamic() {
		double numActivities = 0;
		for (int i = 0; i < edges.size(); i++){
			if (edges.get(i).node1.startsWith("a."))
				numActivities++;
			if (i == edges.size()-1) {
				if (edges.get(i).node2.startsWith("a."))
					numActivities++;
			}
		}
		//this.numActivities = numActivities;
		return numActivities;
	}		
	
	/**
	 * Returns true if this path has a loop, false otherwise
	 * @return true if this path has a loop, false otherwise
	 */
	public boolean hasLoop(){
		boolean hasLoop = false;
		for (int i = 0; i < edges.size(); i++) {
			for (int j = 0; j < edges.size(); j++) {
				if (i != j) {				
					if (edges.get(i).node1 == edges.get(j).node1 || edges.get(i).node1 == edges.get(edges.size()-1).node2) {	
						hasLoop = true;
						break;
					}
				}
			}
		}
		return hasLoop;
	}	
	
	
	@Override	
	public boolean equals(Object object) {

		if (object instanceof Path) {
			Path path = (Path) object;
			
			ArrayList<Edge> pathEdges = path.edges;
			
			if (edges.size()!=pathEdges.size()) {
				return false;
			}
			else {
				for (int i = 0; i < edges.size(); i++) {
					if (!edges.get(i).equals(pathEdges.get(i))) {
						return false;
					}
				}
				return true;
			}
		}
		return false;
	}

	@Override
	public String toString() {
		String s = new String();
		for (int j = 0; j < edges.size(); j++) {	
			s += edges.get(j).node1 + "->";
			if (j == edges.size()-1)
				s += edges.get(j).node2;
		}
		
		//s += "(" + pVal + ")";
		//s += "(" + statistic + ")" + " [" + pVal + "]";
		//s += "(" + numActivities + ")" + " [" + weight + "]";
		//s += "(a=" + numActivities + ", w=" + weight + ", stat=" + statistic + ", pVal=" + pVal + ", ta=" + tempActivities + ")";
		s += "(a=" + numActivities + ", w=" + weight + ", LR=" + round(statistic,2) + ", pVal=" + pVal + ")";
		return s;
	}
	
	public static double round(double value, int places) {
	    if (places < 0) throw new IllegalArgumentException();

	    BigDecimal bd = new BigDecimal(value);
	    bd = bd.setScale(places, BigDecimal.ROUND_HALF_UP);
	    return bd.doubleValue();
	}
	
	/**
	 * Returns whether this path is a subpath of the given path (a path may be a subpath of itself)
	 * @param path
	 * @return whether this path is a subpath of the given path
	 */
	public boolean isSubPathOf (Path path) {
		
		if (this.edges.size() > path.edges.size()) {
			return false;
		}
		else {
			for (int i = 0; i < edges.size(); i++) {
				if (!path.edges.contains(edges.get(i))) {
					return false;
				}
			}		
			return true;
		}
	}
	
	/** 
	 * Returns whether this path contains the given sub path
	 * @param subPath
	 * @return whether this path contains the given sub path
	 */
	public boolean hasSubpath(Path subPath) {	
		return this.edges.containsAll(subPath.edges);
	}

	/**
	 * Returns all sub paths of this path
	 * @return all sub paths of this path
	 */
	public ArrayList<Path> getSubPaths() {
		ArrayList<Path> subPaths = new ArrayList<Path>();
		
		for (int i = 0; i < edges.size(); i++) {
			ArrayList<Edge> pathEdges = new ArrayList<Edge>();
					
			for (int j = i; j < edges.size(); j++) {
				pathEdges.add(edges.get(j));
				
				ArrayList<Edge> temp = new ArrayList<Edge>();
				temp.addAll(pathEdges);
				subPaths.add(new Path(temp));
			}		
		}
		
		return subPaths;
	}
	
	/**
	 * Compares this path to the given path. Returns 1 if this paths has more
	 * activities, -1 if this path has less activities and 0 if both paths
	 * have equal activities.
	 * @param path - the path that this path is being compared to
	 * @return 1 if this paths has more
	 * activities, -1 if this path has less activities and 0 if both paths
	 * have equal activities.
	 */
	public int compareTo(Path path) {
		if (this.numActivities > path.numActivities)
			return 1;
		else if (this.numActivities < path.numActivities)
			return -1;
		else
			return 0;
	}
	
	/**
	 * Returns the node IDs of all nodes in this path
	 * @return the node IDs of all nodes in this path
	 */
	public ArrayList<String> getNodes() {		
		ArrayList<String> nodes = new ArrayList<String>();
		for (int i = 0; i < edges.size(); i++) {
			nodes.add(edges.get(i).node1);
			if (i == edges.size()-1) {
				nodes.add(edges.get(i).node2);
			}
		}
		return nodes;
	}
	
	/**
	 * Returns the activities on this path
	 * @return the activities on this path
	 */
	public ArrayList<Activity> getActivities() {
		ArrayList<Activity> activities = new ArrayList<Activity>();
		for (int i = 0; i < edges.size(); i++) {
			activities.addAll(edges.get(i).activities);			
		}		
		return activities;
	}
	
	
	/**
	 * Returns whether the path contains the given edge
	 * @param edge
	 * @return whether the path contains the given edge
	 */
	public boolean hasEdge(Edge edge) {	
		for (int i = 0; i < edges.size(); i++) {
			if (edges.get(i).equals(edge))
				return true;
		}		
		return false;
	}
	
	
	public static boolean intersects(Path p1, Path p2) {
		//check if p2 contains any of the edges in p1
		for (int i=0; i < p1.edges.size(); i++) {
			if (p2.hasEdge(p1.edges.get(i)))
				return true;
		}
		
		//check if p1 contains any of the edges in p2
		for (int i=0; i < p2.edges.size(); i++) {
			if (p1.hasEdge(p2.edges.get(i)))
				return true;
		}
		
		return false;
	}

	@Override
	protected Path clone() throws CloneNotSupportedException {

		Path path = new Path();
		for (int i = 0; i < this.edges.size(); i++) {
			path.edges.add(this.edges.get(i).clone());
		}
		path.calculatedActivities = this.calculatedActivities;
		path.weight = this.weight;
		path.numActivities = this.numActivities;
		path.clusterActivities = this.clusterActivities;
		path.statistic = this.statistic;
		
		return path;
	}

	
	
}
