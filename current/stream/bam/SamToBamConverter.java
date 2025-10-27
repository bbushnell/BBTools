package stream.bam;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import dna.AminoAcid;
import stream.SamLine;

/**
 * Converts SAM text/SamLine to BAM binary format.
 * Handles encoding of CIGAR, SEQ (4-bit), QUAL, and auxiliary tags.
 *
 * @author Chloe
 * @date October 18, 2025
 */
public class SamToBamConverter {

	private final Map<String, Integer> refMap;
	private final byte[] seqLookup;
	private final Map<Character, Integer> cigarOpMap;

	// SEQ encoding lookup: '=ACMGRSVTWYHKDBN' -> [0,15]
	private static final String SEQ_CHARS = "=ACMGRSVTWYHKDBN";

	// CIGAR operation encoding
	private static final String CIGAR_OPS = "MIDNSHP=X";

	public SamToBamConverter(String[] refNames) {
		// Build reference name to ID map
		refMap = new HashMap<>();
		for (int i = 0; i < refNames.length; i++) {
			refMap.put(refNames[i], i);
		}

		// Build SEQ lookup table
		seqLookup = new byte[256];
		Arrays.fill(seqLookup, (byte)15); // Default to 'N'
		for (int i = 0; i < SEQ_CHARS.length(); i++) {
			char c = SEQ_CHARS.charAt(i);
			seqLookup[c] = (byte)i;
			seqLookup[Character.toLowerCase(c)] = (byte)i;
		}

		// Build CIGAR operation map
		cigarOpMap = new HashMap<>();
		for (int i = 0; i < CIGAR_OPS.length(); i++) {
			cigarOpMap.put(CIGAR_OPS.charAt(i), i);
		}
	}

	/**
	 * Convert a SamLine to BAM binary format.
	 * @return Complete BAM record including block_size prefix
	 */
	public byte[] convertAlignment(SamLine sl) throws IOException {
		ByteArrayOutputStream baos = new ByteArrayOutputStream(64/*sl.estimateBamLength()*/);
		BamWriterHelper writer = new BamWriterHelper(baos);

		// Get reference IDs
		int refID = getRefID(sl.rnameS()); // Use rnameS() which returns String regardless of storage mode
		int nextRefID = getNextRefID(sl.rnext(), refID);

		// Calculate bin
		int bin = calculateBin(sl);

		// Encode CIGAR
		int[] cigarOps = encodeCigar(sl.cigar);

		// Get sequence length
		int seqLen = (sl.seq == null || sl.seq.length == 0) ? 0 : sl.seq.length;

		// Write fixed-length fields (32 bytes)
		writer.writeInt32(refID);
		writer.writeInt32(sl.pos - 1); // Convert 1-based to 0-based
		writer.writeUint8(sl.qname.length() + 1); // Include null terminator
		writer.writeUint8(sl.mapq);
		writer.writeUint16(bin);
		writer.writeUint16(cigarOps.length);
		writer.writeUint16(sl.flag);
		writer.writeUint32(seqLen);
		writer.writeInt32(nextRefID);
		writer.writeInt32(sl.pnext - 1); // Convert 1-based to 0-based
		writer.writeInt32(sl.tlen);

		// Write variable-length fields

		// QNAME with null terminator
		writer.writeString(sl.qname);
		writer.writeUint8(0);

		// CIGAR
		for (int cigOp : cigarOps) {
			writer.writeUint32(cigOp);
		}

		// SEQ (4-bit encoded)
		// SamLine stores reverse-strand sequences in RC form (FLIP_ON_LOAD convention)
		// When writing to BAM, we must RC back to original orientation (matching toBytes() behavior)
		byte[] seqToEncode = sl.seq;
		boolean mapped = refID >= 0;
		boolean reverseStrand = (sl.flag & 0x10) != 0;
		if (mapped && reverseStrand && sl.seq != null && sl.seq.length > 0) {
			seqToEncode = reverseComplement(sl.seq);
		}

		byte[] packedSeq = encodeSeq(seqToEncode);
		writer.writeBytes(packedSeq);

		// QUAL (raw phred scores)
		// Also reverse QUAL for reverse-strand reads (matching toBytes() behavior)
		byte[] qualToEncode = sl.qual;
		if (mapped && reverseStrand && sl.qual != null && sl.qual.length > 0) {
			qualToEncode = reverseArray(sl.qual);
		}
		byte[] rawQual = encodeQual(qualToEncode, seqLen);
		writer.writeBytes(rawQual);

		// Auxiliary tags
		if (sl.optional != null) {
			for (String tag : sl.optional) {
				encodeTag(writer, tag);
			}
		}

		// Return the record bytes without block_size prefix
		// (block_size is added by BamWriter when writing to stream)
		return baos.toByteArray();
	}

	/**
	 * Get reference ID from reference name string.
	 */
	private int getRefID(String rname) {
		if (rname == null || rname.equals("*")) {
			return -1;
		}
		Integer id = refMap.get(rname);
		return (id != null) ? id : -1;
	}

	/**
	 * Get next reference ID from RNEXT.
	 */
	private int getNextRefID(byte[] rnext, int currentRefID) {
		if (rnext == null || rnext.length == 0) {
			return -1;
		}
		if (rnext.length == 1) {
			if (rnext[0] == '*') {
				return -1;
			} else if (rnext[0] == '=') {
				return currentRefID;
			}
		}
		String rnextStr = new String(rnext);
		return getRefID(rnextStr);
	}

	/**
	 * Encode CIGAR string to binary format.
	 * Each operation: (length << 4) | op_code
	 */
	private int[] encodeCigar(String cigar) {
		if (cigar == null || cigar.equals("*")) {
			return new int[0];
		}

		// Parse CIGAR string
		int[] ops = new int[cigar.length()]; // Overestimate
		int count = 0;
		int len = 0;

		for (int i = 0; i < cigar.length(); i++) {
			char c = cigar.charAt(i);
			if (c >= '0' && c <= '9') {
				len = len * 10 + (c - '0');
			} else {
				Integer opCode = cigarOpMap.get(c);
				if (opCode == null) {
					throw new RuntimeException("Unknown CIGAR operation: " + c);
				}
				ops[count++] = (len << 4) | opCode;
				len = 0;
			}
		}

		// Trim to actual size
		return Arrays.copyOf(ops, count);
	}

	/**
	 * Encode SEQ to 4-bit format (2 bases per byte).
	 */
	private byte[] encodeSeq(byte[] seq) {
		if (seq == null || seq.length == 0) {
			return new byte[0];
		}

		byte[] packed = new byte[(seq.length + 1) / 2];
		for (int i = 0; i < seq.length; i++) {
			int val = seqLookup[seq[i] & 0xFF];
			if (i % 2 == 0) {
				packed[i / 2] = (byte)(val << 4);
			} else {
				packed[i / 2] |= (byte)val;
			}
		}
		return packed;
	}

	/**
	 * Encode QUAL to raw phred scores for BAM.
	 * Per SAMv1 spec section 4.2.3: "Base qualities are stored as bytes in the range [0, 93],
	 * without any +33 conversion to printable ASCII."
	 *
	 * NOTE: SamLine already stores qual as raw phred (SamLine.java line 631 does qual[i]-=33),
	 * so we just copy it directly to BAM without further conversion.
	 */
	private byte[] encodeQual(byte[] qual, int seqLen) {
		if (qual == null || qual.length == 0) {
			// Missing quality - use 0xFF per spec
			byte[] missing = new byte[seqLen];
			Arrays.fill(missing, (byte)0xFF);
			return missing;
		}

		// Validate: quality length must match sequence length
		if (qual.length != seqLen) {
			throw new RuntimeException("QUAL length mismatch: qual.length=" + qual.length +
			                           " but seqLen=" + seqLen +
			                           ". SAM data is corrupted.");
		}

		// SamLine.qual is already raw phred - just return a copy
		return qual.clone();
	}

	/**
	 * Encode an auxiliary tag from SAM text format to BAM binary.
	 * Format: TAG:TYPE:VALUE
	 */
	private void encodeTag(BamWriterHelper writer, String tagStr) throws IOException {
		// Parse tag using simple string operations
		if (tagStr.length() < 5) {
			throw new RuntimeException("Invalid tag format: " + tagStr);
		}

		// Extract tag (2 chars), type (1 char), value (rest)
		String tag = tagStr.substring(0, 2);
		char type = tagStr.charAt(3);
		String value = tagStr.substring(5);

		// Write tag (2 bytes)
		writer.writeUint8(tag.charAt(0));
		writer.writeUint8(tag.charAt(1));

		// Write type and value based on type
		switch (type) {
			case 'A': // Printable character
				writer.writeUint8('A');
				writer.writeUint8(value.charAt(0));
				break;

			case 'i': // Integer - choose smallest representation
				long intVal = Long.parseLong(value);
				if (intVal >= Byte.MIN_VALUE && intVal <= Byte.MAX_VALUE) {
					writer.writeUint8('c');
					writer.writeUint8((int)intVal);
				} else if (intVal >= 0 && intVal <= 255) {
					writer.writeUint8('C');
					writer.writeUint8((int)intVal);
				} else if (intVal >= Short.MIN_VALUE && intVal <= Short.MAX_VALUE) {
					writer.writeUint8('s');
					writer.writeInt16((int)intVal);
				} else if (intVal >= 0 && intVal <= 65535) {
					writer.writeUint8('S');
					writer.writeUint16((int)intVal);
				} else if (intVal >= Integer.MIN_VALUE && intVal <= Integer.MAX_VALUE) {
					writer.writeUint8('i');
					writer.writeInt32((int)intVal);
				} else {
					writer.writeUint8('I');
					writer.writeUint32(intVal);
				}
				break;

			case 'f': // Float
				writer.writeUint8('f');
				float floatVal = Float.parseFloat(value);
				writer.writeFloat(floatVal);
				break;

			case 'Z': // String
				writer.writeUint8('Z');
				writer.writeString(value);
				writer.writeUint8(0); // Null terminator
				break;

			case 'H': // Hex string
				writer.writeUint8('H');
				writer.writeString(value);
				writer.writeUint8(0); // Null terminator
				break;

			case 'B': // Array
				writer.writeUint8('B');
				// Parse array: type,val1,val2,...
				String[] parts = value.split(",");
				if (parts.length < 1) {
					throw new RuntimeException("Invalid array tag: " + tagStr);
				}
				char arrayType = parts[0].charAt(0);
				writer.writeUint8(arrayType);
				writer.writeInt32(parts.length - 1); // Count

				for (int i = 1; i < parts.length; i++) {
					switch (arrayType) {
						case 'c':
							writer.writeUint8(Integer.parseInt(parts[i]));
							break;
						case 'C':
							writer.writeUint8(Integer.parseInt(parts[i]));
							break;
						case 's':
							writer.writeInt16(Integer.parseInt(parts[i]));
							break;
						case 'S':
							writer.writeUint16(Integer.parseInt(parts[i]));
							break;
						case 'i':
							writer.writeInt32(Integer.parseInt(parts[i]));
							break;
						case 'I':
							writer.writeUint32(Long.parseLong(parts[i]));
							break;
						case 'f':
							writer.writeFloat(Float.parseFloat(parts[i]));
							break;
						default:
							throw new RuntimeException("Unknown array type: " + arrayType);
					}
				}
				break;

			default:
				throw new RuntimeException("Unknown tag type: " + type);
		}
	}

	/**
	 * Calculate BAM bin using reg2bin algorithm from SAMv1 spec.
	 */
	private int calculateBin(SamLine sl) {
		if (sl.pos <= 0) {
			return 4680; // Special value for unmapped
		}

		int beg = sl.pos - 1; // 0-based
		int end = beg + calculateAlignmentLength(sl.cigar);
		return reg2bin(beg, end);
	}

	/**
	 * Calculate alignment length from CIGAR string.
	 */
	private int calculateAlignmentLength(String cigar) {
		if (cigar == null || cigar.equals("*")) {
			return 0;
		}

		int length = 0;
		int num = 0;

		for (int i = 0; i < cigar.length(); i++) {
			char c = cigar.charAt(i);
			if (c >= '0' && c <= '9') {
				num = num * 10 + (c - '0');
			} else {
				// Operations that consume reference: M, D, N, =, X
				if (c == 'M' || c == 'D' || c == 'N' || c == '=' || c == 'X') {
					length += num;
				}
				num = 0;
			}
		}

		return length;
	}

	/**
	 * Calculate BAM bin for a region [beg, end).
	 * Implementation from SAMv1.pdf page 20.
	 */
	private int reg2bin(int beg, int end) {
		--end;
		if (beg >> 14 == end >> 14) return ((1 << 15) - 1) / 7 + (beg >> 14);
		if (beg >> 17 == end >> 17) return ((1 << 12) - 1) / 7 + (beg >> 17);
		if (beg >> 20 == end >> 20) return ((1 << 9) - 1) / 7 + (beg >> 20);
		if (beg >> 23 == end >> 23) return ((1 << 6) - 1) / 7 + (beg >> 23);
		if (beg >> 26 == end >> 26) return ((1 << 3) - 1) / 7 + (beg >> 26);
		return 0;
	}

	/**
	 * Reverse-complement a DNA sequence.
	 * Used to restore original orientation for reverse-strand reads.
	 */
	private static byte[] reverseComplement(byte[] seq) {
		byte[] rc = new byte[seq.length];
		for (int i = 0, j = seq.length - 1; i < seq.length; i++, j--) {
			rc[i] = AminoAcid.baseToComplementExtended[seq[j]];
		}
		return rc;
	}

	/**
	 * Reverse an array (used for quality scores).
	 */
	private static byte[] reverseArray(byte[] arr) {
		byte[] reversed = new byte[arr.length];
		for (int i = 0, j = arr.length - 1; i < arr.length; i++, j--) {
			reversed[i] = arr[j];
		}
		return reversed;
	}
}
