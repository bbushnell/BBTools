package stream.bam;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.List;

import fileIO.ReadWrite;
import shared.Shared;
import shared.Tools;
import stream.SamLine;
import stream.SamStreamer;
import structures.ListNum;

/**
 * Writes BAM files from SAM data.
 * BAM is the binary compressed version of SAM format.
 *
 * @author Chloe
 * @date October 18, 2025
 */
public class BamWriter {

	private static final boolean DEBUG = false;

	private final OutputStream bgzf;
	private final String filename;
	private final BgzfOutputStream singleThreadedBgzf;
	private final BgzfOutputStreamMT multiThreadedBgzf;
	private BamWriterHelper writer;
	private SamToBamConverter converter;
	private boolean headerWritten = false;
	private boolean closed = false;
	private String[] refNames;
	private int[] refLengths;

	/**
	 * Create a BAM writer for the specified file.
	 */
	public BamWriter(String filename) throws IOException {
		this.filename = filename;
		FileOutputStream fos = new FileOutputStream(filename);
		int zl=Tools.min(ReadWrite.ZIPLEVEL, 6);
		int threads=Tools.mid(1, Shared.threads(), zl>4 ? 16 : zl>3 ? 8 : 4);
		int blockSize = Math.max(1, BgzfSettings.WRITE_BLOCK_SIZE);
		if (BgzfSettings.USE_MULTITHREADED_BGZF && threads>1) {
			multiThreadedBgzf = new BgzfOutputStreamMT(fos, threads, zl, blockSize);
			singleThreadedBgzf = null;
			bgzf = multiThreadedBgzf;
		} else {
			singleThreadedBgzf = new BgzfOutputStream(fos, zl);
			multiThreadedBgzf = null;
			bgzf = singleThreadedBgzf;
		}
		writer = new BamWriterHelper(bgzf);
	}

	/**
	 * Create a BAM writer for the specified file.
	 */
	public BamWriter(OutputStream os) throws IOException {
		this.filename = "os";
		int zl=Tools.min(ReadWrite.ZIPLEVEL, 6);
		int threads=Tools.mid(1, Shared.threads(), zl>4 ? 16 : zl>3 ? 8 : 4);
		int blockSize = Math.max(1, BgzfSettings.WRITE_BLOCK_SIZE);
		if (BgzfSettings.USE_MULTITHREADED_BGZF && threads>1) {
			multiThreadedBgzf = new BgzfOutputStreamMT(os, threads, zl, blockSize);
			singleThreadedBgzf = null;
			bgzf = multiThreadedBgzf;
		} else {
			singleThreadedBgzf = new BgzfOutputStream(os, zl);
			multiThreadedBgzf = null;
			bgzf = singleThreadedBgzf;
		}
		writer = new BamWriterHelper(bgzf);
	}

	/**
	 * Write the BAM header.
	 * Must be called before writeLines().
	 * @param headerLines List of SAM header lines (as byte arrays)
	 */
	public void writeHeader(List<byte[]> headerLines) throws IOException {
		if (headerWritten) {
			throw new RuntimeException("Header already written");
		}

		// Write magic "BAM\1"
		writer.writeBytes(new byte[]{'B', 'A', 'M', 1});

		// Extract references from @SQ lines and build header text
		ArrayList<String> refNamesList = new ArrayList<>();
		ArrayList<Integer> refLengthsList = new ArrayList<>();
		StringBuilder textBuilder = new StringBuilder();

		for (byte[] line : headerLines) {
			// Append to header text
			textBuilder.append(new String(line)).append('\n');

			// Parse @SQ lines to extract reference names and lengths
			if (line.length > 3 && line[0] == '@' && line[1] == 'S' && line[2] == 'Q') {
				String sn = null;
				int ln = 0;

				// Parse using simple string split
				String lineStr = new String(line);
				String[] fields = lineStr.split("\\t");
				for (int i = 1; i < fields.length; i++) {
					if (fields[i].startsWith("SN:")) {
						sn = fields[i].substring(3);
					} else if (fields[i].startsWith("LN:")) {
						ln = Integer.parseInt(fields[i].substring(3));
					}
				}

				if (sn != null && ln > 0) {
					refNamesList.add(sn);
					refLengthsList.add(ln);
				}
			}
		}

		// Convert to arrays
		refNames = refNamesList.toArray(new String[0]);
		refLengths = new int[refLengthsList.size()];
		for (int i = 0; i < refLengths.length; i++) {
			refLengths[i] = refLengthsList.get(i);
		}

		// Write header text section
		byte[] textBytes = textBuilder.toString().getBytes("US-ASCII");
		if (DEBUG) {
			System.err.println("Header: textBytes.length=" + textBytes.length + ", n_ref=" + refNames.length);
		}
		writer.writeUint32(textBytes.length);
		writer.writeBytes(textBytes);

		// Write reference dictionary section
		writer.writeInt32(refNames.length);
		for (int i = 0; i < refNames.length; i++) {
			byte[] nameBytes = refNames[i].getBytes("US-ASCII");
			writer.writeUint32(nameBytes.length + 1); // Include null terminator
			writer.writeBytes(nameBytes);
			writer.writeUint8(0); // Null terminator
			writer.writeUint32(refLengths[i]);
		}

		// Flush the header block
		bgzf.flush();

		// Create converter for alignment records
		converter = new SamToBamConverter(refNames);
		headerWritten = true;
	}

	/**
	 * Write a list of SAM alignments to BAM format.
	 * Can be called multiple times.
	 */
	public void writeLines(List<SamLine> alignments) throws IOException {
		if (!headerWritten) {
			throw new RuntimeException("Must write header first");
		}
		if (closed) {
			throw new IOException("BamWriter already closed");
		}

		for (SamLine sl : alignments) {
			byte[] bamRecord = converter.convertAlignment(sl);
			// Write block_size followed by the record
			writer.writeUint32(bamRecord.length);
			writer.writeBytes(bamRecord);
		}
	}

	/**
	 * Close the BAM file, writing the EOF marker.
	 */
	public void close() throws IOException {
		if (closed) {
			return;
		}
//		System.err.println("Called close");
		if (multiThreadedBgzf != null) {
			multiThreadedBgzf.close();
		} else if (singleThreadedBgzf != null) {
			singleThreadedBgzf.flush();
			singleThreadedBgzf.writeEOF();
			singleThreadedBgzf.close();
		}
		closed = true;
	}

	/**
	 * Generate a BAI index for the written BAM file, using default filename.
	 * Close the writer first to ensure all data is flushed.
	 */
	public void writeIndex() throws IOException {
		writeIndex(null);
	}

	/**
	 * Generate a BAI index for the written BAM file.
	 * @param baiPath Optional explicit index path; defaults to <bam>.bai
	 */
	public void writeIndex(String baiPath) throws IOException {
		if (!closed) {
			throw new IOException("Close BamWriter before writing index");
		}
		if (filename == null) {
			throw new IOException("BamWriter cannot create index without backing filename");
		}
		String target = (baiPath != null && !baiPath.isEmpty()) ? baiPath : filename + ".bai";
		BamIndexWriter.writeIndex(filename, target);
	}

	/**
	 * Main method for testing: convert SAM to BAM.
	 * Usage: java stream.BamWriter input.sam output.bam
	 */
	public static void main(String[] args) throws Exception {
		if (args.length < 2) {
			System.err.println("Usage: java stream.BamWriter input.sam output.bam");
			System.exit(1);
		}

		String samFile = args[0];
		String bamFile = args[1];

		System.err.println("Converting SAM to BAM:");
		System.err.println("  Input:  " + samFile);
		System.err.println("  Output: " + bamFile);

		// Read SAM with SamStreamer
		SamStreamer sls = SamStreamer.makeStreamer(samFile, 4, true, false, -1, false);
		sls.start();

		// Wait for header to be populated
		while (sls.header != null && sls.header.isEmpty()) {
			try {
				Thread.sleep(10);
			} catch (InterruptedException e) {
				break;
			}
		}

		// Make header copy
		List<byte[]> headerCopy = null;
		if (sls.header != null) {
			synchronized (sls.header) {
				headerCopy = new ArrayList<byte[]>(sls.header);
			}
		}

		// Create BAM writer
		BamWriter bw = new BamWriter(bamFile);

		// Write header
		if (headerCopy != null) {
			bw.writeHeader(headerCopy);
		} else {
			throw new RuntimeException("No header found in SAM file");
		}

		// Write alignments
		long count = 0;
		for (ListNum<SamLine> list = sls.nextLines(); list != null; list = sls.nextLines()) {
			bw.writeLines(list.list);
			count += list.size();
		}

		bw.close();

		System.err.println("Converted " + count + " alignments");
		System.err.println("Wrote: " + bamFile);
	}
}
