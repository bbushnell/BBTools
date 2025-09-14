package aligner;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class Parallelogram {
	
	public static void main(String[] args) {
		try {
			convertParallelogramToRectangle(args[0], args[1]);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public static void convertParallelogramToRectangle(String inputFile, String outputFile) throws IOException {
	    // Read input file
	    List<String> lines = new ArrayList<>();
	    BufferedReader reader = new BufferedReader(new FileReader(inputFile));
	    String line;
	    while ((line = reader.readLine()) != null) {
	        if (!line.trim().isEmpty()) {  // Skip empty lines
	            lines.add(line);
	        }
	    }
	    reader.close();
	    
	    if (lines.isEmpty()) {
	        return;  // Nothing to transform
	    }
	    
	    // Calculate dimensions
	    int inputRows = lines.size();
	    int maxWidth = 0;
	    for (String l : lines) {
	        maxWidth = Math.max(maxWidth, l.length());
	    }
	    
	    // The parallelogram should occupy roughly half the matrix
	    // So we'll create a matrix sized for the content
	    char[][] outputMatrix = new char[inputRows][maxWidth];
	    for (char[] row : outputMatrix) {
	        Arrays.fill(row, ' ');
	    }
	    
	    // Process each character in the input
	    for (int i = 0; i < inputRows; i++) {
	        String inputLine = lines.get(i);
	        for (int j = 0; j < inputLine.length(); j++) {
	            char c = inputLine.charAt(j);
	            
	            // Transform coordinates - shift UP by column number
	            int newRow = i - j;
	            
	            // Skip if out of bounds
	            if (newRow < 0 || newRow >= inputRows) {
	                continue;
	            }
	            
	            outputMatrix[newRow][j] = c;
	        }
	    }
	    
	    // Helper method to check if a row is only whitespace
	    boolean[] hasContent = new boolean[inputRows];
	    for (int i = 0; i < inputRows; i++) {
	        for (int j = 0; j < maxWidth; j++) {
	            if (outputMatrix[i][j] != ' ') {
	                hasContent[i] = true;
	                break;
	            }
	        }
	    }
	    
	    // Write output file - only include non-whitespace rows
	    BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile));
	    for (int i = 0; i < inputRows; i++) {
	        if (hasContent[i]) {
	            writer.write(outputMatrix[i]);
	            writer.newLine();
	        }
	    }
	    writer.close();
	}
}
