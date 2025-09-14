package fun;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

public class Crunch {

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
