import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;


public class Utility {

	
	public void parse1 (String allLinesFile, ArrayList<String> subLineFiles) throws IOException {
		
		//BufferedWriter w = new BufferedWriter(new FileWriter(new File("data/output_all_lines.txt")));
		
		
	
		//read in all the file information in arrayLists
		ArrayList<ArrayList<String>> readersData = new ArrayList<ArrayList<String>>();
		for (int i = 0; i < subLineFiles.size(); i++) {
			BufferedReader currentReader = new BufferedReader(new FileReader(new File(subLineFiles.get(i))));
			
			String line = new String();
			ArrayList<String> currentReaderData = new ArrayList<String>();
			while((line = currentReader.readLine()) != null) {
				currentReaderData.add(line);
			}
			readersData.add(currentReaderData);		
		}
		
		//for each line in the master file
		ArrayList<String> output = new ArrayList<String>();
		String line = new String();
		BufferedReader r1 = new BufferedReader(new FileReader(new File(allLinesFile)));
		while((line = r1.readLine()) != null) {
			
			//search for the current line in the other files
			boolean found = false;
			for (int i = 0; i < readersData.size(); i++) {
				
				ArrayList<String> currentFileData = readersData.get(i);
 				for (int j = 0; j < currentFileData.size(); j++) {
 					if (line.equals(currentFileData.get(j))) {
 						found = true;
 						output.add(line + "\t" + (i+1));
 						break;
 					}
 				}		
 				if (found) break;
			}
			if (!found)
				output.add(line + "\t" + -1);
		}
		
		for (int i = 0; i < output.size(); i++) {
			System.out.println(output.get(i));
		}
	}
	
	
	/**
	 * Count the number of times each line in file 1 appears in file 2
	 * @param file1
	 * @param file2
	 * @throws IOException 
	 */
	public void countOccurrences (String file1, String file2) throws IOException {
		
		BufferedReader r1 = new BufferedReader(new FileReader(new File(file1)));
		BufferedReader r2 = new BufferedReader(new FileReader(new File(file2)));
		ArrayList<String> file1Data = new ArrayList<String>();
		ArrayList<String> file2Data = new ArrayList<String>();
		
		//read file 1 into an array
		String line = new String();
		while((line = r1.readLine()) != null) {
			file1Data.add(line);
		}
		
		//read file 2 into an array
		line = new String();
		while((line = r2.readLine()) != null) {
			file2Data.add(line);
		}
		
		for (int i = 0; i < file1Data.size(); i++) {
			int count = 0;
			for (int j = 0; j < file2Data.size(); j++) {
				if (file1Data.get(i).equals(file2Data.get(j)))
					count++;
			}
			System.out.println(file1Data.get(i) + "\t" + count);
		}
		
	}
	
	
	int countNodes (String edgesFile) throws IOException {
		
		ArrayList<String> nodes = new ArrayList<String>();
		
		String line = new String();
		BufferedReader r = new BufferedReader(new FileReader(new File(edgesFile)));
		
		while((line = r.readLine()) != null) {
			
			String[] lineArray = line.split(";");			
			String node1 = lineArray[0].trim();
			String node2 = lineArray[1].trim();
			
			if (!nodes.contains(node1))
				nodes.add(node1);
			if (!nodes.contains(node2))
				nodes.add(node2);			
		}
		
		return nodes.size();
	}
	
	
	public static void main (String[] args) throws IOException {
		
		Utility u = new Utility();
		System.out.println("node count: " + u.countNodes("C:\\data\\pedestrian_fatalities_oc\\experimental_data\\edges.txt"));
		
		System.out.println("done");
	}
}
