package stream.bam;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;

/**
 * Converts BAM binary alignment records to SAM text format.
 * Handles all BAM field encodings including 4-bit SEQ, packed CIGAR, and auxiliary tags.
 *
 * @author Chloe
 * @date October 18, 2025
 */
public class BamToSamConverter {

	private final String[] refNames;

	// Lookup tables for decoding
	private static final char[] SEQ_LOOKUP = "=ACMGRSVTWYHKDBN".toCharArray();
	private static final char[] CIGAR_OPS = "MIDNSHP=X".toCharArray();

	public BamToSamConverter(String[] refNames) {
		this.refNames = refNames;
	}

	/**
	 * Convert a BAM alignment record to SAM text format.
	 * @param bamRecord The raw BAM record bytes (not including block_size field)
	 * @return Tab-delimited SAM line as byte array
	 */
	public byte[] convertAlignment(byte[] bamRecord) {
		ByteBuffer bb = ByteBuffer.wrap(bamRecord).order(ByteOrder.LITTLE_ENDIAN);

		// Read fixed-length fields
		int refID = bb.getInt();
		int pos = bb.getInt();
		int l_read_name = bb.get() & 0xFF;
		int mapq = bb.get() & 0xFF;
		int bin = bb.getShort() & 0xFFFF; // Ignore bin
		int n_cigar_op = bb.getShort() & 0xFFFF;
		int flag = bb.getShort() & 0xFFFF;
		long l_seq = bb.getInt() & 0xFFFFFFFFL;
		int next_refID = bb.getInt();
		int next_pos = bb.getInt();
		int tlen = bb.getInt();

		// Read variable-length fields
		byte[] readNameBytes = new byte[l_read_name];
		bb.get(readNameBytes);
		String qname = new String(readNameBytes, 0, l_read_name - 1); // Exclude NUL terminator

		// Decode CIGAR
		StringBuilder cigar = new StringBuilder();
		if (n_cigar_op == 0) {
			cigar.append('*');
		} else {
			for (int i = 0; i < n_cigar_op; i++) {
				int cigOp = bb.getInt();
				int opLen = cigOp >>> 4;
				int op = cigOp & 0xF;
				cigar.append(opLen).append(CIGAR_OPS[op]);
			}
		}

		// Decode SEQ (4-bit encoded, 2 bases per byte)
		String seq;
		if (l_seq == 0) {
			seq = "*";
		} else {
			StringBuilder seqBuilder = new StringBuilder((int)l_seq);
			int numBytes = (int)((l_seq + 1) / 2);
			for (int i = 0; i < numBytes; i++) {
				int packed = bb.get() & 0xFF;
				seqBuilder.append(SEQ_LOOKUP[packed >>> 4]);
				if (i * 2 + 1 < l_seq) {
					seqBuilder.append(SEQ_LOOKUP[packed & 0xF]);
				}
			}
			seq = seqBuilder.toString();
		}

		// Decode QUAL (raw phred scores, add 33 for SAM)
		String qual;
		if (l_seq == 0) {
			qual = "*";
		} else {
			byte[] qualBytes = new byte[(int)l_seq];
			bb.get(qualBytes);

			// Check if QUAL is missing (all 0xFF)
			boolean allFF = true;
			for (int i = 0; i < qualBytes.length; i++) {
				if (qualBytes[i] != (byte)0xFF) {
					allFF = false;
					break;
				}
			}

			if (allFF) {
				qual = "*";
			} else {
				// Convert to ASCII (phred + 33)
				for (int i = 0; i < qualBytes.length; i++) {
					qualBytes[i] = (byte)(qualBytes[i] + 33);
				}
				try {
					qual = new String(qualBytes, 0, qualBytes.length, "US-ASCII");
				} catch (java.io.UnsupportedEncodingException e) {
					throw new RuntimeException("US-ASCII not supported", e); // Should never happen
				}
			}
		}

		// Decode auxiliary tags
		StringBuilder tags = new StringBuilder();
		while (bb.hasRemaining()) {
			byte[] tagBytes = new byte[2];
			bb.get(tagBytes);
			String tag = new String(tagBytes, 0, 2);
			char type = (char)(bb.get() & 0xFF);

			if (tags.length() > 0) {
				tags.append('\t');
			}
			tags.append(tag).append(':');

			switch (type) {
				case 'A': // Printable character
					tags.append("A:").append((char)(bb.get() & 0xFF));
					break;
				case 'c': // int8_t
					tags.append("i:").append((int)bb.get());
					break;
				case 'C': // uint8_t
					tags.append("i:").append(bb.get() & 0xFF);
					break;
				case 's': // int16_t
					tags.append("i:").append((int)bb.getShort());
					break;
				case 'S': // uint16_t
					tags.append("i:").append(bb.getShort() & 0xFFFF);
					break;
				case 'i': // int32_t
					tags.append("i:").append(bb.getInt());
					break;
				case 'I': // uint32_t
					long uintVal = bb.getInt() & 0xFFFFFFFFL;
					tags.append("i:").append(uintVal);
					break;
				case 'f': // float
					tags.append("f:").append(bb.getFloat());
					break;
				case 'Z': // Null-terminated string
					StringBuilder str = new StringBuilder();
					byte b;
					while ((b = bb.get()) != 0) {
						str.append((char)(b & 0xFF));
					}
					tags.append("Z:").append(str);
					break;
				case 'H': // Hex string
					StringBuilder hex = new StringBuilder();
					while ((b = bb.get()) != 0) {
						hex.append((char)(b & 0xFF));
					}
					tags.append("H:").append(hex);
					break;
				case 'B': // Array
					char arrayType = (char)(bb.get() & 0xFF);
					int count = bb.getInt();
					tags.append("B:").append(arrayType);
					for (int i = 0; i < count; i++) {
						tags.append(',');
						switch (arrayType) {
							case 'c':
								tags.append((int)bb.get());
								break;
							case 'C':
								tags.append(bb.get() & 0xFF);
								break;
							case 's':
								tags.append((int)bb.getShort());
								break;
							case 'S':
								tags.append(bb.getShort() & 0xFFFF);
								break;
							case 'i':
								tags.append(bb.getInt());
								break;
							case 'I':
								tags.append(bb.getInt() & 0xFFFFFFFFL);
								break;
							case 'f':
								tags.append(bb.getFloat());
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

		// Build SAM line
		// QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL [TAGS]
		StringBuilder sam = new StringBuilder();

		// QNAME
		sam.append(qname).append('\t');

		// FLAG
		sam.append(flag).append('\t');

		// RNAME
		if (refID < 0 || refID >= refNames.length) {
			sam.append('*');
		} else {
			sam.append(refNames[refID]);
		}
		sam.append('\t');

		// POS (BAM is 0-based, SAM is 1-based)
		sam.append(pos + 1).append('\t');

		// MAPQ
		sam.append(mapq).append('\t');

		// CIGAR
		sam.append(cigar).append('\t');

		// RNEXT
		if (next_refID < 0) {
			sam.append('*');
		} else if (next_refID == refID) {
			sam.append('=');
		} else if (next_refID < refNames.length) {
			sam.append(refNames[next_refID]);
		} else {
			sam.append('*');
		}
		sam.append('\t');

		// PNEXT (BAM is 0-based, SAM is 1-based)
		sam.append(next_pos + 1).append('\t');

		// TLEN
		sam.append(tlen).append('\t');

		// SEQ
		sam.append(seq).append('\t');

		// QUAL
		sam.append(qual);

		// Tags (if any)
		if (tags.length() > 0) {
			sam.append('\t').append(tags);
		}

		// SAM format requires US-ASCII encoding (SAM spec page 1)
		try {
			return sam.toString().getBytes("US-ASCII");
		} catch (java.io.UnsupportedEncodingException e) {
			throw new RuntimeException("US-ASCII not supported", e); // Should never happen
		}
	}

	/**
	 * Reverse-complement a DNA sequence string.
	 * Used to restore original sequence orientation for reverse-strand reads.
	 */
	private static String reverseComplement(String seq) {
		StringBuilder rc = new StringBuilder(seq.length());
		for (int i = seq.length() - 1; i >= 0; i--) {
			char base = seq.charAt(i);
			char comp;
			switch (base) {
				case 'A': comp = 'T'; break;
				case 'T': comp = 'A'; break;
				case 'G': comp = 'C'; break;
				case 'C': comp = 'G'; break;
				case 'a': comp = 't'; break;
				case 't': comp = 'a'; break;
				case 'g': comp = 'c'; break;
				case 'c': comp = 'g'; break;
				case 'N': comp = 'N'; break;
				case 'n': comp = 'n'; break;
				default: comp = base; // Keep ambiguity codes as-is
			}
			rc.append(comp);
		}
		return rc.toString();
	}
}
