package stream.bam;

import java.io.*;
import java.util.*;

import shared.LineParser1;
import stream.BamLineStreamer;
import stream.SamLine;

/**
 * Minimal test case for BAM round-trip bug.
 * Tests writing and reading back a minimal SAM alignment.
 *
 * @author Chloe
 * @date October 18, 2025
 */
public class BamRoundTripTest {

	public static void main(String[] args) throws Exception {
		System.out.println("=== BAM Round-Trip Bug Test ===\n");

		// Create minimal SAM header
		List<byte[]> header = new ArrayList<>();
		header.add("@HD\tVN:1.4\tSO:unsorted".getBytes());
		header.add("@SQ\tSN:chr1\tLN:1000".getBytes());

		// Create minimal SAM alignment (unmapped read - simplest case)
		String samLine = "read1\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGTAC\t??????????";
		System.out.println("Original SAM line:");
		System.out.println(samLine);
		System.out.println("Tabs: " + countTabs(samLine) + " (expected: 10)");
		System.out.println();

		// Parse to SamLine
		SamLine sl = new SamLine(new LineParser1('\t').set(samLine.getBytes()));

		// Write to BAM
		String bamFile = "test_minimal.bam";
		BamWriter writer = new BamWriter(bamFile);
		writer.writeHeader(header);

		// Get the BAM record bytes for debugging
		SamToBamConverter converter = new SamToBamConverter(new String[]{"chr1"});
		byte[] bamRecord = converter.convertAlignment(sl);

		System.out.println("BAM record created:");
		System.out.println("Total length: " + bamRecord.length + " bytes");
		System.out.println("First 64 bytes (hex):");
		hexDump(bamRecord, 0, Math.min(64, bamRecord.length));
		System.out.println();

		// Write the alignment
		writer.writeLines(Arrays.asList(sl));
		writer.close();

		System.out.println("Wrote BAM file: " + bamFile);
		System.out.println();

		// Read back with BamLineStreamer
		System.out.println("Reading back from BAM...");
		try {
			BamLineStreamer bls = new BamLineStreamer(bamFile, 1, true, true, -1, false);
			bls.start();

			structures.ListNum<SamLine> list = bls.nextLines();
			if (list != null && list.size() > 0) {
				SamLine readBack = list.get(0);

				// Convert back to text for comparison
				byte[] samBytes = readBack.toText().toBytes();
				String reconverted = new String(samBytes);

				System.out.println("SUCCESS! Read back SAM line:");
				System.out.println(reconverted);
				System.out.println("Tabs: " + countTabs(reconverted) + " (expected: 10)");
				System.out.println();

				// Analyze the intermediate SAM text from BamToSamConverter
				System.out.println("Analyzing BAM→SAM conversion:");
				BamToSamConverter bamToSam = new BamToSamConverter(new String[]{"chr1"});

				// Extract just the alignment record (skip block_size prefix)
				byte[] recordOnly = Arrays.copyOfRange(bamRecord, 4, bamRecord.length);
				byte[] samFromBam = bamToSam.convertAlignment(recordOnly);
				String samText = new String(samFromBam);

				System.out.println("SAM from BAM converter:");
				System.out.println(samText);
				System.out.println("Length: " + samFromBam.length + " bytes");
				System.out.println("Tabs: " + countTabs(samText));
				System.out.println();

				// Character-by-character analysis
				System.out.println("Character analysis:");
				for (int i = 0; i < Math.min(samFromBam.length, 200); i++) {
					byte b = samFromBam[i];
					if (b == '\t') {
						System.out.print("[TAB]");
					} else if (b == '\n') {
						System.out.print("[LF]");
					} else if (b >= 32 && b < 127) {
						System.out.print((char)b);
					} else {
						System.out.print("[" + (b & 0xFF) + "]");
					}
				}
				System.out.println();

			} else {
				System.out.println("ERROR: No alignments read back!");
			}

//			bls.close();

		} catch (Exception e) {
			System.out.println("\nERROR during read-back:");
			e.printStackTrace();

			// Hex dump the BAM file for analysis
			System.out.println("\nBAM file hex dump (first 512 bytes):");
			FileInputStream fis = new FileInputStream(bamFile);
			byte[] fileBytes = new byte[512];
			int n = fis.read(fileBytes);
			fis.close();
			hexDump(fileBytes, 0, n);
		}
	}

	private static int countTabs(String s) {
		int count = 0;
		for (int i = 0; i < s.length(); i++) {
			if (s.charAt(i) == '\t') count++;
		}
		return count;
	}

	private static void hexDump(byte[] data, int offset, int length) {
		for (int i = offset; i < offset + length; i++) {
			if ((i - offset) % 16 == 0) {
				if (i > offset) System.out.println();
				System.out.printf("%04X: ", i - offset);
			}
			System.out.printf("%02X ", data[i] & 0xFF);
		}
		System.out.println();
	}
}
