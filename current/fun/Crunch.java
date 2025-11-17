package fun;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

/**
 * File size aggregation utility for processing file metadata from a
 * pipe-delimited input file. Calculates total size of files matching
 * specific criteria and converts the result to terabytes.
 * Filters files containing 'F' in the 7th column of the metadata.
 *
 * @author Brian Bushnell
 */
public class Crunch {

	/**
	 * Program entry point that processes a pipe-delimited file to calculate
	 * total size of files matching filtering criteria.
	 * Reads file metadata from pipe-delimited text file, filters files
	 * containing 'F' in the 7th column, aggregates total file size across
	 * matching files, and converts final size to terabytes.
	 *
	 * @param args Command-line arguments where args[0] is the input file path
	 */
	public static void main(String[] args) {
		try {
			Scanner scanner = new Scanner(new File(args[0]));
			long total_size = 0;

			while (scanner.hasNextLine()) {
				String line = scanner.nextLine();
				String[] splitLine = line.split("\\|");
				if(!splitLine[6].contains("F")) {
					continue;
				}
				long size = Long.parseLong(splitLine[3]);
				long atime = Long.parseLong(splitLine[11]);
				long mtime = Long.parseLong(splitLine[12]);

				total_size+=size;
			}

			System.out.printf("Total Size: %d TB = " +total_size, (total_size / 1024 / 1024 / 1024 / 1024));
			scanner.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}

}
