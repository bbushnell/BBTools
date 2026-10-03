package prot;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.nio.file.attribute.BasicFileAttributes;
import java.nio.file.attribute.FileTime;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;

/**
 * Constants and shared helpers for the Increment 3B member-index artifact
 * ({@code member_index.bin}, design v1-v7, root-accepted
 * {@code root_acceptance_design_v7.md}). Standalone: introduces no dependency on, and does
 * not modify, the already-accepted Increment 1/2 {@code HbmBundleFormat}/
 * {@code HbmBundleBuilder}/{@code HbmBundleLoader}/{@code HbmMemberManifest} classes --
 * per the acceptance's "do not apply shared code" scope, this file and its siblings are a
 * self-contained new artifact family, not an extension of the existing one.
 *
 * <p>All multi-byte integers are big-endian with no alignment padding, matching this
 * project's established convention. The frozen ID hash (design v3 sec2, v4 sec2): the first
 * 8 bytes of SHA-256(id_bytes), interpreted as a big-endian {@code i64}.</p>
 *
 * @author Eru
 */
public final class HbmMemberIndexFormat {

	private HbmMemberIndexFormat(){}

	/*--------------------------------------------------------------*/
	/*----------------      Magic and version       ----------------*/
	/*--------------------------------------------------------------*/

	private static final byte[] MAGIC={(byte)'M', (byte)'Q', (byte)'M', (byte)'I'};//"MQMI"
	public static final int MAGIC_LEN=4;
	public static byte[] magic(){return MAGIC.clone();}
	public static boolean magicMatches(final byte[] data, final int off){
		checkRange(data, off, MAGIC_LEN, "magicMatches");
		for(int i=0; i<MAGIC_LEN; i++){if(data[off+i]!=MAGIC[i]){return false;}}
		return true;
	}
	public static final int FORMAT_VERSION=1;

	/*--------------------------------------------------------------*/
	/*----------------      Fixed-header layout      ----------------*/
	/*--------------------------------------------------------------*/
	// magic[4] | format_version(u16) | entry_count(i64) | capacity(i64) | load_factor_permille(u16) | arena_chunk_bytes(i32)
	public static final int OFF_MAGIC=0;
	public static final int OFF_FORMAT_VERSION=4;
	public static final int OFF_ENTRY_COUNT=6;
	public static final int OFF_CAPACITY=14;
	public static final int OFF_LOAD_FACTOR_PERMILLE=22;
	public static final int OFF_ARENA_CHUNK_BYTES=24;
	public static final int FIXED_HEADER_LEN=28;

	public static final int SHA256_LEN=32;
	public static final int LOAD_FACTOR_PERMILLE_DEFAULT=700;//0.70

	/** Distinct from the Increment 1/2 bundle's {@code CAP_REPID} (which governs family
	 *  representative IDs only) -- this cap governs arbitrary source/member IDs (design v3
	 *  sec2, v4 sec2). */
	public static final int CAP_SOURCE_ID_BYTES=2048;

	/** Reviewed ceiling on a single arena chunk (slice-1 v2, root_review_slice1_v1.md sec3):
	 *  128 MiB -- generous for any real ID, bounded so a malformed {@code arena_chunk_bytes}
	 *  field cannot trigger an unbounded allocation before any record is even read. */
	public static final int MAX_ARENA_CHUNK_BYTES=128<<20;

	/** Fixed per-slot payload of the four parallel primitive arrays in
	 *  {@link HbmMemberIndexTable} (long hash + int rank + long arenaLoc + int idLen). */
	public static final int TABLE_BYTES_PER_SLOT=8+4+8+4;

	public static final String PASS1_FORMAT_ID="hbm_3b_pass1_stage_v1";

	/** Frozen buffer size used by {@link #streamingWholeFileSha256} -- recorded canonically in
	 *  the Pass-1 manifest (slice-1 v3, root_review_slice1_v3.md sec5) so the resource record
	 *  names the real memory cost of hashing, not just the table/arena budgets. */
	public static final int HASH_BUFFER_BYTES=1<<16;

	/*--------------------------------------------------------------*/
	/*----------------          Helpers              ----------------*/
	/*--------------------------------------------------------------*/

	private static void checkRange(final byte[] data, final int off, final int len, final String who){
		if(data==null){throw new IllegalArgumentException(who+": null data");}
		if(off<0 || len<0 || off>data.length || len>data.length-off){
			throw new IllegalArgumentException(who+": bad range off="+off+" len="+len+" size="+data.length);
		}
	}

	public static byte[] sha256(final byte[] data, final int off, final int len){
		checkRange(data, off, len, "sha256");
		final MessageDigest md;
		try{md=MessageDigest.getInstance("SHA-256");}
		catch(NoSuchAlgorithmException e){throw new RuntimeException("SHA-256 unavailable", e);}
		md.update(data, off, len);
		return md.digest();
	}

	/** Frozen hash: first 8 bytes of SHA-256(id_bytes), big-endian i64 (design v3 sec2). */
	public static long idHash(final byte[] idBytes, final int off, final int len){
		final byte[] d=sha256(idBytes, off, len);
		long h=0;
		for(int i=0; i<8; i++){h=(h<<8)|(d[i]&0xffL);}
		return h;
	}

	/** The one production hasher every real build/load must use. Tests inject a different
	 *  {@link HbmIdHasher} (e.g. a constant hash) against the same table code to prove
	 *  correctness does not depend on this choice. */
	public static final HbmIdHasher FROZEN_HASHER=new HbmIdHasher(){
		@Override
		public long hash(byte[] idBytes, int off, int len){return idHash(idBytes, off, len);}
	};

	/** Smallest power of two &gt;= x, for x&gt;0. Overflow-checked. */
	public static long nextPowerOfTwo(final long x){
		if(x<=0){throw new IllegalArgumentException("nextPowerOfTwo requires x>0, got "+x);}
		if(x>(1L<<62)){throw new IllegalArgumentException("nextPowerOfTwo overflow risk for x="+x);}
		long p=1;
		while(p<x){p<<=1;}
		return p;
	}

	/** Overflow-checked total primitive-array payload bytes for a table of this capacity
	 *  (slice-1 v2, root_review_slice1_v1.md sec3: computed before allocation, checked
	 *  against a caller-supplied reviewed budget). */
	public static long tableBytesFor(final long capacity){
		if(capacity<=0){throw new IllegalArgumentException("capacity must be positive: "+capacity);}
		if(capacity>Long.MAX_VALUE/TABLE_BYTES_PER_SLOT){throw new IllegalArgumentException("capacity too large for overflow-checked sizing: "+capacity);}
		return capacity*TABLE_BYTES_PER_SLOT;
	}

	/** Table capacity from an exact tracked-row count at a 0.70 load factor (design v3 sec1/
	 *  v4 sec1: sized from Pass 0.5's MEASURED count, never an estimate), overflow-checked. */
	public static long capacityFor(final long exactTrackedRowCount){
		if(exactTrackedRowCount<0){throw new IllegalArgumentException("negative row count: "+exactTrackedRowCount);}
		if(exactTrackedRowCount==0){throw new IllegalArgumentException("cannot size a table for zero tracked rows");}
		if(exactTrackedRowCount>Long.MAX_VALUE/10){throw new IllegalArgumentException("row count too large for overflow-checked sizing: "+exactTrackedRowCount);}
		//capacity >= count/0.70 = count*10/7, ceiling division.
		final long minCapacity=(exactTrackedRowCount*10+6)/7;
		return nextPowerOfTwo(minCapacity);
	}

	/** Streaming whole-file SHA-256 (slice-1 v3, root_review_slice1_v2.md sec4): reads the file
	 *  through a bounded buffer via {@code FileChannel}, never {@code Files.readAllBytes} --
	 *  the whole point being that a large sealed artifact (e.g. {@code member_index.bin}) is
	 *  never materialized in memory solely to hash it. Returns lowercase hex directly since
	 *  every call site compares against a canonical hex manifest field. */
	public static String streamingWholeFileSha256(final String path) throws IOException{
		final MessageDigest md;
		try{md=MessageDigest.getInstance("SHA-256");}
		catch(NoSuchAlgorithmException e){throw new RuntimeException("SHA-256 unavailable", e);}
		try(FileChannel ch=FileChannel.open(Paths.get(path), StandardOpenOption.READ)){
			final long size=ch.size();
			final ByteBuffer buf=ByteBuffer.allocate(HASH_BUFFER_BYTES);
			long pos=0;
			while(pos<size){
				buf.clear();
				final long rem=size-pos;
				if(rem<buf.capacity()){buf.limit((int)rem);}
				final int r=ch.read(buf, pos);
				if(r<0){throw new IOException("EOF hashing whole file at "+pos+" of "+size+": "+path);}
				if(r==0){throw new IOException("NO_PROGRESS_READ: FileChannel.read returned 0 at "+pos+": "+path);}
				buf.flip();
				md.update(buf);
				pos+=r;
			}
			return toHexLower(md.digest());
		}
	}

	/** A capturable snapshot of an input file's identity (slice-1 v3/v4,
	 *  root_review_slice1_v3.md and root_review_slice1_v4.md item 1: "do not treat two
	 *  successive opens of the same pathname as proof of the same bytes" -- and path+byte-count
	 *  ALONE is not enough either; root proved a same-size content rewrite passes that check).
	 *  Canonical path (symlinks resolved), byte count, the filesystem's file key (device+inode
	 *  on POSIX, {@code null} if the platform/filesystem does not expose one), and the file's
	 *  last-modified time -- captured immediately before a scan and required unchanged
	 *  immediately after (via {@link #requireUnchangedIdentity}). This metadata layer is
	 *  necessary but production call sites additionally recompute and compare the accepted
	 *  whole-file hash after the parsing scan, before any reusable seal is written (see
	 *  {@code HbmPass05Stage.write()}/{@code HbmMemberIndexBuilder}) -- the hash recheck is the
	 *  definitive, content-based closure; this identity snapshot is a cheap, immediate fast-fail
	 *  layer in front of it, and the mechanism this class alone can prove (mtime changing on
	 *  rewrite) independent of any hash recomputation. */
	public static final class InputIdentity{
		public final String canonicalPath;
		public final long byteCount;
		public final Object fileKey;//may be null: not every platform/filesystem exposes one
		public final FileTime lastModifiedTime;
		InputIdentity(final String canonicalPath, final long byteCount, final Object fileKey, final FileTime lastModifiedTime){
			this.canonicalPath=canonicalPath; this.byteCount=byteCount;
			this.fileKey=fileKey; this.lastModifiedTime=lastModifiedTime;
		}
	}

	public static InputIdentity captureInputIdentity(final String path) throws IOException{
		final Path real=Paths.get(path).toRealPath();
		final BasicFileAttributes attrs=Files.readAttributes(real, BasicFileAttributes.class);
		return new InputIdentity(real.toString(), attrs.size(), attrs.fileKey(), attrs.lastModifiedTime());
	}

	/** Requires {@code path}'s CURRENT identity to exactly match a previously captured one --
	 *  call immediately before and immediately after each scan of the same file. The file-key
	 *  comparison is skipped when either snapshot lacks one (platform does not support it, not
	 *  a mismatch); path, byte count, and modification time are always compared. */
	public static void requireUnchangedIdentity(final InputIdentity before, final String path, final String tag) throws IOException{
		final InputIdentity after=captureInputIdentity(path);
		if(!before.canonicalPath.equals(after.canonicalPath)){
			throw new RuntimeException(tag+": input path identity changed across scans (before="+before.canonicalPath+
				", after="+after.canonicalPath+")");
		}
		if(before.byteCount!=after.byteCount){
			throw new RuntimeException(tag+": input identity changed across scans (before="+before.byteCount+
				" bytes, after="+after.byteCount+" bytes)");
		}
		if(before.fileKey!=null && after.fileKey!=null && !before.fileKey.equals(after.fileKey)){
			throw new RuntimeException(tag+": input file key changed across scans (before="+before.fileKey+
				", after="+after.fileKey+") -- the path now names a different underlying file");
		}
		if(!before.lastModifiedTime.equals(after.lastModifiedTime)){
			throw new RuntimeException(tag+": input identity changed across scans (before mtime="+before.lastModifiedTime+
				", after mtime="+after.lastModifiedTime+")");
		}
	}

	public static String toHexLower(final byte[] digest){
		if(digest==null || digest.length!=SHA256_LEN){
			throw new IllegalArgumentException("expected 32-byte digest, got "+(digest==null?"null":String.valueOf(digest.length)));
		}
		final char[] hex="0123456789abcdef".toCharArray();
		final char[] out=new char[digest.length*2];
		for(int i=0; i<digest.length; i++){
			final int b=digest[i]&0xff;
			out[i*2]=hex[b>>>4];
			out[i*2+1]=hex[b&0xf];
		}
		return new String(out);
	}
}
