package fun;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 * Analyzes disk space usage by recursively scanning directories and
 * generating an HTML report showing directory sizes sorted by size.
 * Useful for identifying large directories consuming disk space.
 * @author Brian Bushnell
 */
public class DiskSpaceAnalyzer {

	/** Internal class to track directory size information.
	 * Contains path and total bytes for a directory. */
	static class DirSize {
		String path;
		long bytes;

		DirSize(String path, long bytes) {
			this.path = path;
			this.bytes = bytes;
		}
	}

	private static Comparator<DirSize> SIZE_COMPARATOR = new Comparator<DirSize>() {
		@Override
		public int compare(DirSize a, DirSize b) {
			return Long.compare(b.bytes, a.bytes);
		}
	};

	public static void analyzeDiskSpace(String rootPath) throws IOException {
		File root = new File(rootPath);
		List<DirSize> dirSizes = new ArrayList<>();

		analyzeDirRecursive(root, dirSizes);

		// Sort by size, descending
		Collections.sort(dirSizes, new Comparator<DirSize>() {
			@Override
			public int compare(DirSize a, DirSize b) {
				return Long.compare(b.bytes, a.bytes);
			}
		});

		generateHTMLReport(dirSizes, rootPath);
	}

	private static long analyzeDirRecursive(File dir, List<DirSize> results) {
		long total = 0;
		if (!dir.exists() || !dir.isDirectory()) return total;

		File[] files = dir.listFiles();
		if (files == null) return total;

		for (File file : files) {
			if (file.isDirectory()) {
				long dirSize = analyzeDirRecursive(file, results);
				total += dirSize;
				results.add(new DirSize(file.getPath(), dirSize));
			} else {
				total += file.length();
			}
		}

		return total;
	}

	private static void generateHTMLReport(List<DirSize> dirSizes, String rootPath) throws IOException {
		DecimalFormat df = new DecimalFormat("#.##");

		try (FileWriter writer = new FileWriter("disk_usage.html")) {
			writer.write("<!DOCTYPE html>\n<html>\n<head>\n");
			writer.write("<title>Disk Usage Report</title>\n");
			writer.write("<style>body { font-family: Arial; }</style>\n");
			writer.write("</head>\n<body>\n");
			writer.write("<h1>Disk Usage Report for " + rootPath + "</h1>\n");
			writer.write("<table border='1'>\n");
			writer.write("<tr><th>Path</th><th>Size (MB)</th><th>Percentage</th></tr>\n");

			Collections.sort(dirSizes, SIZE_COMPARATOR);
			long totalSize=0;
			for (DirSize dir : dirSizes) {totalSize+=dir.bytes;}

			for (DirSize dir : dirSizes) {
				double sizeInMB = dir.bytes / (1024.0 * 1024.0);
				double percentage = (dir.bytes * 100.0) / totalSize;

				writer.write(String.format(
						"<tr><td>%s</td><td>%.2f</td><td>%.2f%%</td></tr>\n", 
						dir.path, sizeInMB, percentage
						));
			}

			writer.write("</table>\n</body>\n</html>");
		}

		System.out.println("Disk usage report generated: disk_usage.html");
	}

	public static void main(String[] args) throws IOException {
		String path = args.length > 0 ? args[0] : ".";
		analyzeDiskSpace(path);
	}
}
