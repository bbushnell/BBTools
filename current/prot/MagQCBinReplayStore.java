package prot;

import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Random;

import parse.LineParser1;
import structures.ByteBuilder;
import structures.IntList;
import structures.LongList;

/**
 * Disk-backed, ordinal-indexed selections for one model in a shared manifest.
 * Every input row is validated before any replay output opens. Only selected
 * payloads reach the private spool; heap retains primitive offsets rather than
 * contig lists. Canonical pool_N IDs have constant-size uniqueness bookkeeping;
 * irregular IDs require an exact fallback set during preflight.
 *
 * <p>The private binary spool contains a four-byte length followed by one TSV
 * row. Source text uses the canonical ByteFile reader. JDK buffered binary
 * output and random access provide the seekable scratch format; row parsing
 * and formatting remain the manifest's existing implementations.</p>
 * @author Yoimiya
 */
public final class MagQCBinReplayStore implements AutoCloseable {

	/** Consumer-specific source and row capability checks run for the entire input. */
	public interface Validator {
		void header(MagQCBinManifest.Manifest header);
		void bin(MagQCBinManifest.Bin bin,MagQCBinManifest.Manifest header);
	}

	/** Creates an owned temporary spool; null directory uses java.io.tmpdir. */
	public static MagQCBinReplayStore build(String file,String model,Path directory,Validator validator){
		return build(file,model,directory,validator,null);
	}

	/** Creates an owned temporary spool, optionally applying a reproducible per-split shuffle. */
	public static MagQCBinReplayStore build(String file,String model,Path directory,Validator validator,Long shuffleSeed){
		return build(file,model,directory,validator,shuffleSeed,false);
	}

	/** Creates a replay store with optional in-memory raw rows for shuffle-heavy workers. */
	public static MagQCBinReplayStore build(String file,String model,Path directory,Validator validator,
			Long shuffleSeed,boolean memoryRows){
		MagQCBinManifest.validateModelName(model);
		if(validator==null){throw new IllegalArgumentException("Replay requires full-manifest validation");}
		try{return buildChecked(file,model,directory,validator,shuffleSeed,memoryRows);}
		catch(IOException e){throw new RuntimeException("Could not build replay spool for "+file,e);}
	}

	/** Validates, writes selected payloads and constructs contiguous per-split offset arrays. */
	private static MagQCBinReplayStore buildChecked(String file,String model,Path directory,Validator validator,Long shuffleSeed,
			boolean memoryRows) throws IOException{
		Path spool=null;
		try{
			final MagQCBinManifest.Manifest header;
			final IntList[] ordinals={new IntList(),new IntList()};
			final LongList[] positions={new LongList(),new LongList()};
			final ArrayList<Retained> retainedTrain=memoryRows ? new ArrayList<Retained>() : null;
			final ArrayList<Retained> retainedVal=memoryRows ? new ArrayList<Retained>() : null;
			final Ids ids=new Ids();
			long bytes=0,rows=0;
			int maxLength=0;
			try(MagQCBinManifest.Reader reader=new MagQCBinManifest.Reader(file)){
				header=reader.header(); validator.header(header);
				spool=(directory==null ? Files.createTempFile("magqc-replay-",".bin") :
					Files.createTempFile(directory,"magqc-replay-",".bin"));
				try(DataOutputStream output=new DataOutputStream(new BufferedOutputStream(Files.newOutputStream(spool),65536))){
					final ByteBuilder row=new ByteBuilder(4096);
					for(MagQCBinManifest.Bin bin=reader.next(); bin!=null; bin=reader.next()){
						ids.add(bin.id); rows++;
						validator.bin(bin,header);
						final int ordinal=MagQCBinManifest.modelOrdinal(bin,model);
						if(ordinal<0){continue;}
						final int split=bin.split.equals("train") ? 0 : 1;
						ordinals[split].add(ordinal); positions[split].add(bytes);
						row.clear(); MagQCBinManifest.append(bin,row);
						final int length=row.length()-1; // The binary record does not need the terminal newline.
						assert(length>0) : "A formatted manifest row contains eleven nonempty fields";
						if(memoryRows){(split==0 ? retainedTrain : retainedVal).add(new Retained(ordinal,Arrays.copyOf(row.array,length)));}
						output.writeInt(length); output.write(row.array,0,length);
						bytes=Math.addExact(bytes,4L+length); maxLength=Math.max(maxLength,length);
					}
				}
			}
			final Order[] ordered={order(ordinals[0],positions[0],model,"train",shuffleSeed,0),
				order(ordinals[1],positions[1],model,"val",shuffleSeed,1)};
			if(ordered[0].offsets.length==0 && ordered[1].offsets.length==0){throw new IllegalArgumentException("No manifest rows for binmodel="+model);}
			final byte[][][] memory=memoryRows ? retainedRows(retainedTrain,retainedVal,ordered) : null;
			return new MagQCBinReplayStore(spool,header,model,
				new long[][]{ordered[0].offsets,ordered[1].offsets},
				new int[][]{ordered[0].sourceOrdinals,ordered[1].sourceOrdinals},memory,maxLength,rows,bytes);
		}catch(IOException | RuntimeException | Error failure){
			if(spool!=null){try{Files.deleteIfExists(spool);}catch(IOException cleanup){failure.addSuppressed(cleanup);}}
			throw failure;
		}
	}

	private static byte[][][] retainedRows(ArrayList<Retained> retainedTrain,ArrayList<Retained> retainedVal,Order[] ordered){
		final ArrayList<Retained>[] retained=new ArrayList[]{retainedTrain,retainedVal};
		final byte[][][] result=new byte[2][][];
		for(int split=0; split<2; split++){
			final byte[][] byOrdinal=new byte[retained[split].size()][];
			for(Retained row : retained[split]){
				if(row.ordinal<0 || row.ordinal>=byOrdinal.length || byOrdinal[row.ordinal]!=null){
					throw new IllegalArgumentException("Retained replay row ordinal is not dense: "+row.ordinal);
				}
				byOrdinal[row.ordinal]=row.bytes;
			}
			result[split]=new byte[ordered[split].sourceOrdinals.length][];
			for(int i=0; i<result[split].length; i++){result[split][i]=byOrdinal[ordered[split].sourceOrdinals[i]];}
		}
		return result;
	}

	private static final class Retained {
		Retained(int ordinal_,byte[] bytes_){ordinal=ordinal_; bytes=bytes_;}
		final int ordinal;
		final byte[] bytes;
	}

	/** Rejects holes and duplicate ordinals without allocating to an untrusted maximum ordinal. */
	private static Order order(IntList ordinals,LongList positions,String model,String split,Long shuffleSeed,int splitIndex){
		assert(ordinals.size==positions.size) : "Every selected ordinal records exactly one spool position";
		final long[] orderedOffsets=new long[ordinals.size];
		final int[] orderedOrdinals=new int[ordinals.size];
		final long[] byOrdinal=new long[ordinals.size];
		Arrays.fill(byOrdinal,-1);
		for(int i=0; i<byOrdinal.length; i++){
			final int ordinal=ordinals.get(i);
			if(ordinal<0 || ordinal>=byOrdinal.length || byOrdinal[ordinal]>=0){
				throw new IllegalArgumentException("Noncontiguous or duplicate row "+ordinal+" for model "+model+" split "+split);
			}
			byOrdinal[ordinal]=positions.get(i);
		}
		for(int i=0; i<orderedOffsets.length; i++){
			orderedOffsets[i]=byOrdinal[i]; orderedOrdinals[i]=i;
		}
		if(shuffleSeed!=null && orderedOffsets.length>1){
			final Random random=new Random(shuffleSeed.longValue() ^ (splitIndex==0 ? 0x6a09e667f3bcc909L : 0xbb67ae8584caa73bL));
			for(int i=orderedOffsets.length-1; i>0; i--){
				final int j=random.nextInt(i+1);
				final long offset=orderedOffsets[i]; orderedOffsets[i]=orderedOffsets[j]; orderedOffsets[j]=offset;
				final int ordinal=orderedOrdinals[i]; orderedOrdinals[i]=orderedOrdinals[j]; orderedOrdinals[j]=ordinal;
			}
		}
		return new Order(orderedOffsets,orderedOrdinals);
	}

	private static final class Order {
		Order(long[] offsets_,int[] sourceOrdinals_){offsets=offsets_; sourceOrdinals=sourceOrdinals_;}
		final long[] offsets;
		final int[] sourceOrdinals;
	}

	private MagQCBinReplayStore(Path file,MagQCBinManifest.Manifest header,String model,long[][] offsets,
			int[][] sourceOrdinals_,byte[][][] memoryRows_,int maxLength,long rows,long bytes) throws IOException{
		this.file=file; this.header=header; this.model=model; this.offsets=offsets;
		sourceOrdinals=sourceOrdinals_;
		memoryRows=memoryRows_;
		this.maxLength=maxLength; rowsRead=rows; spoolBytes=bytes;
		input=(memoryRows==null ? new RandomAccessFile(file.toFile(),"r") : null);
	}

	/** Empty validated metadata, retaining the declared realization split policy. */
	public MagQCBinManifest.Manifest header(){return header;}
	/** Number of selected rows; order is model ordinal order unless a shuffle seed was supplied. */
	public int size(String split){return offsets[splitIndex(split)].length;}
	/** Number of all source rows validated, including rows belonging only to other models. */
	public long rowsRead(){return rowsRead;}
	/** Bytes of private selected-bin payload retained on disk rather than heap. */
	public long spoolBytes(){return spoolBytes;}

	/** Reads only the requested row and checks its split/membership before returning it. */
	public MagQCBinManifest.Bin get(String split,int ordinal){
		if(closed){throw new IllegalStateException("Replay spool is closed");}
		final long[] index=offsets[splitIndex(split)];
		final int[] source=sourceOrdinals[splitIndex(split)];
		if(ordinal<0 || ordinal>=index.length){throw new IndexOutOfBoundsException("Replay row "+ordinal+" of "+index.length);}
		try{
			final byte[] bytes;
			if(memoryRows!=null){bytes=memoryRows[splitIndex(split)][ordinal];}
			else{
				input.seek(index[ordinal]);
				final int length=input.readInt();
				if(length<1 || length>maxLength){throw new IOException("Invalid replay record length: "+length);}
				bytes=new byte[length]; input.readFully(bytes);
			}
			parser.set(bytes);
			final MagQCBinManifest.Bin bin=MagQCBinManifest.parseBin(parser,header.schemaVersion,file.toString());
			if(!bin.split.equals(split) || MagQCBinManifest.modelOrdinal(bin,model)!=source[ordinal]){
				throw new IOException("Replay offset points to the wrong model row: "+ordinal);
			}
			return bin;
		}catch(IOException e){throw new RuntimeException("Could not read replay spool "+file,e);}
	}

	/** Opens an independent seek/parser cursor for one replay worker. */
	Reader openReader(){return new Reader(file,header,model,offsets,sourceOrdinals,memoryRows,maxLength);}

	/** Thread-confined random-access replay cursor; metadata and offsets are shared read-only. */
	static final class Reader implements AutoCloseable {
		private final Path file;
		private final MagQCBinManifest.Manifest header;
		private final String model;
		private final long[][] offsets;
		private final int[][] sourceOrdinals;
		private final byte[][][] memoryRows;
		private final int maxLength;
		private final RandomAccessFile input;
		private final LineParser1 parser=new LineParser1((byte)'\t');
		private boolean closed;

		private Reader(Path file_,MagQCBinManifest.Manifest header_,String model_,long[][] offsets_,int[][] sourceOrdinals_,
				byte[][][] memoryRows_,int maxLength_){
			file=file_; header=header_; model=model_; offsets=offsets_; sourceOrdinals=sourceOrdinals_; memoryRows=memoryRows_; maxLength=maxLength_;
			try{input=(memoryRows==null ? new RandomAccessFile(file.toFile(),"r") : null);}
			catch(IOException e){throw new RuntimeException("Could not open replay worker cursor "+file,e);}
		}

		MagQCBinManifest.Bin get(String split,int ordinal){
			if(closed){throw new IllegalStateException("Replay worker cursor is closed");}
			final long[] index=offsets[splitIndex(split)];
			final int[] source=sourceOrdinals[splitIndex(split)];
			if(ordinal<0 || ordinal>=index.length){throw new IndexOutOfBoundsException("Replay row "+ordinal+" of "+index.length);}
			try{
				final byte[] bytes;
				if(memoryRows!=null){bytes=memoryRows[splitIndex(split)][ordinal];}
				else{
					input.seek(index[ordinal]);
					final int length=input.readInt();
					if(length<1 || length>maxLength){throw new IOException("Invalid replay record length: "+length);}
					bytes=new byte[length]; input.readFully(bytes);
				}
				parser.set(bytes);
				final MagQCBinManifest.Bin bin=MagQCBinManifest.parseBin(parser,header.schemaVersion,file.toString());
				if(!bin.split.equals(split) || MagQCBinManifest.modelOrdinal(bin,model)!=source[ordinal]){
					throw new IOException("Replay offset points to the wrong model row: "+ordinal);
				}
				return bin;
			}catch(IOException e){throw new RuntimeException("Could not read replay spool "+file,e);}
		}

		@Override public void close(){
			if(closed){return;}
			closed=true;
			try{if(input!=null){input.close();}}catch(IOException e){throw new RuntimeException("Could not close replay worker cursor "+file,e);}
		}
	}

	/** Closes and deletes only this store's private temporary file. */
	@Override public void close(){
		if(closed){return;}
		closed=true;
		IOException failure=null;
		try{if(input!=null){input.close();}}catch(IOException e){failure=e;}
		try{Files.delete(file);}catch(IOException e){if(failure==null){failure=e;}else{failure.addSuppressed(e);}}
		if(failure!=null){throw new RuntimeException("Could not close/delete replay spool "+file,failure);}
	}

	/** Refuses unsupported splits rather than silently treating them as validation. */
	private static int splitIndex(String split){
		if("train".equals(split)){return 0;}
		if("val".equals(split)){return 1;}
		throw new IllegalArgumentException("Expected train or val: "+split);
	}

	/** Exact uniqueness: a canonical seen prefix needs no individual String entries. */
	static final class Ids {
		void add(String id){
			final long number=poolNumber(id);
			if((number>=0 && number<nextPoolId) || irregular.contains(id)){
				throw new IllegalArgumentException("Duplicate bin_id "+id);
			}
			if(number==nextPoolId){nextPoolId++;}
			else{irregular.add(id);}
		}
		/** Only the merger's canonical decimal spelling qualifies for prefix compression. */
		private static long poolNumber(String id){
			if(!id.startsWith("pool_") || id.length()==5 || (id.length()>6 && id.charAt(5)=='0')){return -1;}
			long value=0;
			for(int i=5; i<id.length(); i++){
				final int digit=id.charAt(i)-'0';
				if(digit<0 || digit>9 || value>(Long.MAX_VALUE-digit)/10){return -1;}
				value=value*10+digit;
			}
			return value;
		}
		long nextPoolId;
		final HashSet<String> irregular=new HashSet<String>();
	}

	private final Path file;
	private final MagQCBinManifest.Manifest header;
	private final String model;
	private final long[][] offsets;
	private final int[][] sourceOrdinals;
	private final byte[][][] memoryRows;
	private final int maxLength;
	private final long rowsRead,spoolBytes;
	private final RandomAccessFile input;
	private final LineParser1 parser=new LineParser1((byte)'\t');
	private boolean closed;
}
