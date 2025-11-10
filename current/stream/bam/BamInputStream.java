package stream.bam;

import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Objects;

import fileIO.FileFormat;
import shared.Shared;
import stream.SamLine;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * Presents BAM content as SAM text through the {@link InputStream} interface.
 * The implementation wraps {@link Streamer} so callers can treat a BAM
 * file as a plain-text stream of SAM header lines followed by SAM alignment
 * lines, one per newline.
 *
 * <p>This is intentionally inefficient (it round-trips through {@link SamLine})
 * but convenient for components that expect a text stream.</p>
 */
public class BamInputStream extends InputStream {

	public BamInputStream(String fname, boolean ordered) {
		this(FileFormat.testInput(fname, FileFormat.BAM, null, true, false), ordered, 
			Math.min(6, Shared.threads()/2+1));
	}

	public BamInputStream(String fname, boolean ordered, int threads) {
		this(FileFormat.testInput(fname, FileFormat.BAM, null, true, false), ordered, threads);
	}

	public BamInputStream(FileFormat ff, boolean ordered, int threads) {
		Objects.requireNonNull(ff, "FileFormat must not be null");
		streamer = StreamerFactory.makeSamOrBamStreamer(ff, threads, true, ordered, -1L, false);
		streamer.start();
		headerQueue = new ArrayDeque<>();
	}

	@Override
	public int read() throws IOException {
		if (!ensureBuffer()) {
			return -1;
		}
		return currentBuffer[currentOffset++] & 0xFF;
	}

	@Override
	public int read(byte[] b, int off, int len) throws IOException {
		if (b == null) {
			throw new NullPointerException();
		}
		if (off < 0 || len < 0 || len > b.length - off) {
			throw new IndexOutOfBoundsException();
		}
		if (len == 0) {
			return 0;
		}
		if (!ensureBuffer()) {
			return -1;
		}

		int copied = 0;
		while (copied < len) {
			if (currentOffset >= currentBuffer.length) {
				if (!ensureBuffer()) {
					break;
				}
			}
			int toCopy = Math.min(len - copied, currentBuffer.length - currentOffset);
			System.arraycopy(currentBuffer, currentOffset, b, off + copied, toCopy);
			currentOffset += toCopy;
			copied += toCopy;
		}

		return copied == 0 ? -1 : copied;
	}

	@Override
	public synchronized void close() {
		if(closed) {return;}
		streamer.close();
		closed = true;
		currentList = null;
		currentBuffer = null;
	}

	private boolean ensureBuffer() {
		if (closed) {
			return false;
		}
		if (currentBuffer != null && currentOffset < currentBuffer.length) {
			return true;
		}
		currentBuffer = null;
		currentOffset = 0;

		// Header lines first
		if (!headerLoaded) {
			ArrayList<byte[]> sharedHeader = stream.SamReadInputStream.getSharedHeader(false);
			if (sharedHeader == null) {
				for (int attempts = 0; attempts < 100 && sharedHeader == null; attempts++) {
					try {
						Thread.sleep(10L);
					} catch (InterruptedException ignored) {
						Thread.currentThread().interrupt();
					}
					sharedHeader = stream.SamReadInputStream.getSharedHeader(false);
				}
				if (sharedHeader == null) {
					sharedHeader = stream.SamReadInputStream.getSharedHeader(true);
				}
			}
			if (sharedHeader != null) {
				headerQueue.addAll(sharedHeader);
			}
			headerLoaded = true;
		}
		if (!headerQueue.isEmpty()) {
			byte[] headerLine = headerQueue.poll();
			if (headerLine != null) {
				currentBuffer = appendNewline(headerLine);
				return true;
			}
		}

		while (true) {
			if (currentList == null || samIndex >= currentList.size()) {
				currentList = streamer.nextLines();
				samIndex = 0;
				if (currentList == null) {
					closed = true;
					return false;
				}
				if (currentList.size() == 0) {
					continue;
				}
			}

			ArrayList<SamLine> samLines = currentList.list;
			if (samIndex < samLines.size()) {
				SamLine sam = samLines.get(samIndex++);
				if (sam == null) {
					continue;
				}
				ByteBuilder bb = sam.toText();
				bb.append('\n');
				currentBuffer = bb.toBytes();
				return true;
			}
		}
	}

	private static byte[] appendNewline(byte[] line) {
		int len = line.length;
		while (len > 0 && (line[len - 1] == '\r' || line[len - 1] == '\n')) {
			len--;
		}
		ByteBuilder bb = new ByteBuilder(len + 1);
		bb.append(line, 0, len);
		bb.append('\n');
		return bb.toBytes();
	}

	private final Streamer streamer;
	private final ArrayDeque<byte[]> headerQueue;
	private boolean headerLoaded = false;

	private ListNum<SamLine> currentList = null;
	private int samIndex = 0;

	private byte[] currentBuffer = null;
	private int currentOffset = 0;

	private boolean closed = false;
}
