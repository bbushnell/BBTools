package stream.bam;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;

import shared.LineParser1;
import stream.SamLine;
import structures.ByteBuilder;

/**
 * {@link OutputStream} facade that accepts SAM text and emits BAM output.
 * Each newline-terminated SAM line written to this stream is parsed into a
 * {@link SamLine} and handed to {@link SamToBamConverter}, so callers can
 * treat BAM as a drop-in text sink.
 */
public class BamOutputStream extends OutputStream {

	public BamOutputStream(String filename) throws IOException {
		this(new FileOutputStream(filename), true);
	}
	
	public BamOutputStream(String filename, int compression, int threads) throws IOException {
		this(new FileOutputStream(filename), true, compression, threads);
	}

	public BamOutputStream(OutputStream out, boolean closeUnderlying) {
		this(out, closeUnderlying, BgzfSettings.WRITE_COMPRESSION_LEVEL, BgzfSettings.WRITE_THREADS);
	}

	public BamOutputStream(OutputStream out, boolean closeUnderlying, int compression, int threads) {
		this.closeUnderlying = closeUnderlying;
		if (BgzfSettings.USE_MULTITHREADED_BGZF && threads>1) {
			int blockSize = Math.max(1, BgzfSettings.WRITE_BLOCK_SIZE);
			mtOut = new BgzfOutputStreamMT(out, threads, compression, blockSize);
			stOut = null;
			bgzf = mtOut;
//			System.err.println("MT: compression="+compression+", threads="+threads+", blockSize="+blockSize);
		} else {
			stOut = new BgzfOutputStream(out, compression);
			mtOut = null;
			bgzf = stOut;
//			System.err.println("ST: compression="+compression+", threads="+threads);
		}
		helper = new BamWriterHelper(bgzf);
	}

	@Override
	public void write(int b) throws IOException {
		ensureOpen();
		if (b == '\n') {
			flushLine();
		} else {
			lineBuffer.append((byte) b);
		}
	}

	@Override
	public void write(byte[] b, int off, int len) throws IOException {
		ensureOpen();
		int end = off + len;
		for (int i = off; i < end; i++) {
			byte ch = b[i];
			if (ch == '\n') {
				flushLine();
			} else if (ch != '\r') {
				lineBuffer.append(ch);
			}
		}
	}

	@Override
	public void flush() throws IOException {
		ensureOpen();
		flushLine();
		flushSamBuffer();
		bgzf.flush();
	}

	@Override
	public void close() throws IOException {
		if (closed) {
			return;
		}
		try {
			flushLine();
			flushSamBuffer();
			if (!headerWritten) {
				writeHeader();
			}
			if (mtOut != null) {
				mtOut.close();
			} else if (stOut != null) {
				stOut.writeEOF();
				stOut.close();
			} else if (closeUnderlying) {
				bgzf.close();
			} else {
				bgzf.flush();
			}
		} finally {
			closed = true;
		}
	}

	private void flushLine() throws IOException {
		if (lineBuffer.isEmpty()) {
			return;
		}
		byte[] raw = lineBuffer.toBytes();
		lineBuffer.clear();
		int len = raw.length;
		while (len > 0 && (raw[len - 1] == '\r' || raw[len - 1] == '\n')) {
			len--;
		}
		if (len == 0) {
			return;
		}
		byte[] line = (len == raw.length) ? raw : trimBytes(raw, len);
		if (!headerWritten && line[0] == '@') {
			headerLines.add(line);
			return;
		}
		if (!headerWritten) {
			writeHeader();
		}
		SamLine sam = new SamLine(new LineParser1('\t').set(line));
		samBuffer.add(sam);
		if (samBuffer.size() >= FLUSH_THRESHOLD) {
			flushSamBuffer();
		}
	}

	private void flushSamBuffer() throws IOException {
		if (!headerWritten || samBuffer.isEmpty()) {
			return;
		}
		for (SamLine sam : samBuffer) {
			byte[] record = converter.convertAlignment(sam);
			helper.writeUint32(record.length);
			helper.writeBytes(record);
		}
		samBuffer.clear();
	}

	private void writeHeader() throws IOException {
		if (headerWritten) {
			return;
		}
		helper.writeBytes(new byte[]{'B', 'A', 'M', 1});

		StringBuilder textBuilder = new StringBuilder();
		List<String> refs = new ArrayList<>();
		List<Integer> refLengths = new ArrayList<>();
		for (byte[] line : headerLines) {
			String lineStr = new String(line, StandardCharsets.US_ASCII);
			textBuilder.append(lineStr).append('\n');
			if (lineStr.startsWith("@SQ")) {
				String[] fields = lineStr.split("\\t");
				String sn = null;
				Integer ln = null;
				for (int i = 1; i < fields.length; i++) {
					String field = fields[i];
					if (field.startsWith("SN:")) {
						sn = field.substring(3);
					} else if (field.startsWith("LN:")) {
						ln = Integer.parseInt(field.substring(3));
					}
				}
				if (sn != null && ln != null) {
					refs.add(sn);
					refLengths.add(ln);
				}
			}
		}

		byte[] headerText = textBuilder.toString().getBytes(StandardCharsets.US_ASCII);
		helper.writeUint32(headerText.length);
		helper.writeBytes(headerText);

		helper.writeInt32(refs.size());
		for (int i = 0; i < refs.size(); i++) {
			byte[] nameBytes = refs.get(i).getBytes(StandardCharsets.US_ASCII);
			helper.writeUint32(nameBytes.length + 1);
			helper.writeBytes(nameBytes);
			helper.writeUint8(0);
			helper.writeUint32(refLengths.get(i));
		}

		if (converter == null) {
			converter = new SamToBamConverter(refs.toArray(new String[0]));
		}

		headerWritten = true;
	}

	private void ensureOpen() throws IOException {
		if (closed) {
			throw new IOException("Stream closed");
		}
	}

	private static byte[] trimBytes(byte[] data, int len) {
		byte[] trimmed = new byte[len];
		System.arraycopy(data, 0, trimmed, 0, len);
		return trimmed;
	}

	private static final int FLUSH_THRESHOLD = 256;

	private final boolean closeUnderlying;
	private final BgzfOutputStreamMT mtOut;
	private final BgzfOutputStream stOut;
	private final OutputStream bgzf;
	private final BamWriterHelper helper;

	private final ArrayList<byte[]> headerLines = new ArrayList<>();
	private final ArrayList<SamLine> samBuffer = new ArrayList<>();
	private final ByteBuilder lineBuffer = new ByteBuilder(4096);

	private SamToBamConverter converter;
	private boolean headerWritten = false;
	private boolean closed = false;
}
