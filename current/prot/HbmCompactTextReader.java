package prot;

import java.io.IOException;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.ArrayDeque;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import parse.LineParser1;
import shared.Shared;
import structures.ByteBuilder;

/** Reads self-contained compact HBM counts and optional per-model annotations.
 * All rows and headers are covered by one LF-normalized sha80 footer. The
 * entire file is validated before graphs are published to scoring callers.
 * @author Collei
 */
public final class HbmCompactTextReader {

	public static final String FORMAT="hbm_compact_v1";
	/** Diagnostic switch; normal loading uses bounded per-family byte blocks. */
	static boolean BLOCK_INPUT=true;
	static final String CUTOFF_UNITS="identity_percent\tblosum_raw\thbm_path_relative";
	static boolean matches(Path path){
		final String name=path.toString().toLowerCase(Locale.ROOT);
		return name.endsWith(".hbmc") || name.endsWith(".hbmc.gz");
	}

	/** A standalone reader needs no separate consensus, family list or provenance file. */
	public static Result read(String path) throws IOException{
		return read(path, Math.min(256, Shared.threads()), null);
	}

	/** Bounded family tasks preserve file order; no model escapes before all checks pass. */
	static Result read(String path, int threads, Timings timings) throws IOException{
		if(threads<1 || threads>256){throw bad("threads must be in [1,256]");}
		final long start=timings==null ? 0 : System.nanoTime();
		final Input input=new Input(path, BLOCK_INPUT, timings); boolean complete=false;
		final ArrayList<String> ids=new ArrayList<String>();
		final ArrayList<AAGraph> graphs=new ArrayList<AAGraph>();
		final ArrayList<HbmCompactMetadata> metadata=new ArrayList<HbmCompactMetadata>();
		final HashSet<String> names=new HashSet<String>(), headers=new HashSet<String>();
		final String[] provenance=new String[16];
		final ExecutorService pool=threads==1 ? null : Executors.newFixedThreadPool(threads);
		final ArrayDeque<Future<Family>> pending=new ArrayDeque<Future<Family>>();
		int families=-1; String filter=null;
		try{
			long reading=timings==null ? 0 : System.nanoTime();
			long ioBefore=0, hashBefore=0;
			byte[] line=input.next();
			require(line, "#format\t"+FORMAT);
			line=input.next();
			for(; line!=null && !starts(line, "#name\t"); line=input.next()){
				input.fields.set(line); final LineParser1 row=input.fields;
				final String key=row.parseString(0);
				if(!headers.add(key)){throw bad("duplicate global header "+key);}
				if(key.equals("#graph_contract")){require(line, key+"\t"+HbmDenseTextBundle.CONTRACT);}
				else if(key.equals("#cutoff_units")){require(line, key+"\t"+CUTOFF_UNITS);}
				else if(key.equals("#coordinates")){require(line, key+"\tzero_based_inclusive");}
				else if(key.equals("#sums")){require(line, key+"\toriginal");}
				else if(key.equals("#filter")){
					filter=field(row);
					if(!(filter.equals("none") || filter.equals("singletons") || filter.equals("states") || filter.equals("both"))){throw bad("unknown filter");}
				}else if(key.equals("#n_families")){
					field(row); families=row.parseInt(1);
					if(families<1 || families>100000){throw bad("family count outside supported bounds");}
				}else if(key.startsWith("#provenance_")){
					final int index=Integer.parseInt(key.substring(12));
					if(index<0 || index>=16 || !key.equals("#provenance_"+index)){throw bad("invalid provenance index");}
					provenance[index]=DigestSuffix.requireSuffix(field(row), key);
				}else if(key.equals("#columns")){
					if(row.terms()!=5){throw bad("invalid column description");}
				}else{throw bad("unknown global header "+key);}
			}
			for(String key:new String[]{"#graph_contract", "#cutoff_units", "#coordinates", "#sums", "#filter", "#n_families", "#columns"}){
				if(!headers.contains(key)){throw bad("missing "+key);}
			}
			for(String pin:provenance){if(pin==null){throw bad("missing semantic provenance");}}
			final boolean lossy=filter.equals("singletons") || filter.equals("both");
			long nodes=0;
			if(timings!=null){
				final long elapsed=System.nanoTime()-reading;
				timings.readNanos+=elapsed; timings.headerNanos+=elapsed-input.lineNanos-input.hashNanos;
			}
			for(int rank=0; rank<families; rank++){
				reading=timings==null ? 0 : System.nanoTime();
				ioBefore=input.lineNanos; hashBefore=input.hashNanos;
				input.fields.set(required(line)); final LineParser1 row=input.fields;
				if(!row.termEquals("#name", 0)){throw bad("missing family name at rank "+rank);}
				final String encoded=field(row), id=IdentifierCodec.decode(encoded);
				if(!IdentifierCodec.encode(id).equals(encoded) || id.isEmpty() || !names.add(id)){throw bad("invalid or duplicate family name");}
				row.set(required(input.next()));
				if(!row.termEquals("#alphabet", 0)){throw bad("missing alphabet");}
				final int[] alphabet=alphabet(field(row));
				row.set(required(input.next()));
				if(!row.termEquals("#consensus", 0) || row.terms()!=2){throw bad("missing consensus");}
				final byte[] pivot=Blosum62.encode(row.parseByteArray(1), id);
				if(pivot.length<1 || pivot.length>HbmBundleFormat.CAP_L){throw bad("consensus length outside supported bounds");}
				final HbmCompactMetadata info=new HbmCompactMetadata();
				line=required(input.next()); row.set(line);
				while(!row.termEquals("#rows", 0)){
					info.parse(row); line=required(input.next()); row.set(line);
				}
				field(row); final int count=row.parseInt(1);
				if(count<pivot.length || count>HbmBundleFormat.CAP_TOTAL_NODES ||
						(nodes+=count+(long)pivot.length)>HbmBundleFormat.CAP_TOTAL_NODES){throw bad("invalid or excessive node count");}
				if(info.start>=pivot.length || info.stop>=pivot.length){throw bad("endpoint bound outside consensus");}
				if(timings!=null){timings.headerNanos+=System.nanoTime()-reading-(input.lineNanos-ioBefore)-(input.hashNanos-hashBefore);}
				final Rows lines=input.rows(count);
				ids.add(id); metadata.add(info);
				line=input.next();
				if(timings!=null){timings.readNanos+=System.nanoTime()-reading;}
				if(pool==null){append(family(lines, alphabet, pivot, lossy, timings!=null), graphs, timings);}
				else{
					pending.addLast(pool.submit(() -> family(lines, alphabet, pivot, lossy, timings!=null)));
					if(pending.size()>=2*threads){append(get(pending.removeFirst()), graphs, timings);}
				}
			}
			input.fields.set(required(line));
			if(!input.fields.termEquals("#checksum", 0) || !field(input.fields).equals(DigestSuffix.fromDigest(input.digest.digest()))){throw bad("missing or incorrect checksum");}
			if(input.next()!=null){throw bad("trailing data");}
			while(!pending.isEmpty()){append(get(pending.removeFirst()), graphs, timings);}
			complete=true;
			return new Result(ids.toArray(new String[0]), graphs.toArray(new AAGraph[0]),
				metadata.toArray(new HbmCompactMetadata[0]), provenance, filter);
		}finally{
			if(pool!=null){
				pool.shutdownNow(); boolean interrupted=false;
				for(;;){
					try{if(pool.awaitTermination(1, TimeUnit.DAYS)){break;}}
					catch(InterruptedException e){interrupted=true;}
				}
				if(interrupted){Thread.currentThread().interrupt();}
			}
			if(timings!=null){timings.totalNanos=System.nanoTime()-start;}
			if(timings!=null){
				timings.lineNanos=input.lineNanos; timings.hashNanos=input.hashNanos;
				timings.rawIoNanos=input.raw==null ? 0 : input.raw.ioNanos;
			}
			if(input.close() && complete){throw new IOException("Compact input close failed");}
		}
	}

	/** One primitive buffer separates numeric decoding from graph allocation without per-cell objects. */
	private static Family family(Rows rows, int[] alphabet, byte[] pivot, boolean lossy, boolean timed){
		final long splitStart=timed ? System.nanoTime() : 0;
		final ArrayList<byte[]> lines=rows.lines();
		final long start=timed ? System.nanoTime() : 0;
		final int stride=3+alphabet.length;
		if((long)lines.size()*stride>Integer.MAX_VALUE-8){throw bad("decoded family exceeds array cap");}
		final int[] decoded=new int[lines.size()*stride];
		final LineParser1 row=new LineParser1('\t');
		int anchor=-1, chain=0;
		for(int r=0, pos=0; r<lines.size(); r++, pos+=stride){
			if(Thread.currentThread().isInterrupted()){throw bad("family parsing interrupted");}
			row.set(lines.get(r));
			final boolean insertion=row.termEquals('+', 0); final int offset=insertion ? 1 : 0;
			if(insertion){
				if(anchor<0 || anchor+1>=pivot.length || ++chain>HbmBundleFormat.CAP_CHAIN){throw bad("invalid insertion anchor/chain");}
			}else{
				if(++anchor>=pivot.length){throw bad("too many reference rows");}
				chain=0;
			}
			decoded[pos]=insertion ? 1 : 0;
			final int sum=HbmDenseTextLoader.integer(row, offset); decoded[pos+1]=sum;
			int end=row.terms();
			if(end>offset+1 && row.termStartsWith("-", end-1)){
				if(insertion || anchor==0 || anchor==pivot.length-1){throw bad("invalid deletion position");}
				final int n=HbmDenseTextLoader.integer(row, --end, 1);
				if(n<1){throw bad("zero deletion must be omitted");}
				decoded[pos+2]=n;
			}
			final int width=end-offset-1;
			if(sum<1 || width<0 || width>alphabet.length){throw bad("invalid row sum/count width");}
			long actual=0;
			for(int c=0; c<width; c++){
				final int n=HbmDenseTextLoader.integer(row, offset+1+c);
				decoded[pos+3+c]=n; actual+=n;
			}
			if(actual>sum || (!lossy && actual!=sum)){throw bad("residue count sum disagrees with declared filtering");}
		}
		if(anchor!=pivot.length-1){throw bad("missing reference rows");}
		final long parsed=timed ? System.nanoTime() : 0;
		final AAGraph graph=new AAGraph(pivot, 0); HbmBundleLoader.checkGraphKnobs(graph);
		anchor=-1; AAGraphNode tail=null;
		for(int pos=0; pos<decoded.length; pos+=stride){
			if(Thread.currentThread().isInterrupted()){throw bad("family reconstruction interrupted");}
			final AAGraphNode node;
			if(decoded[pos]!=0){
				assert(tail!=null) : "Validated compact INS rows must follow a REF anchor";
				node=new AAGraphNode(Blosum62.X_CODE, AAGraphNode.INS, anchor+1); tail.insEdge=node;
			}else{
				node=graph.ref[++anchor]; Arrays.fill(node.count, 0); Arrays.fill(node.weight, 0);
				graph.del[anchor].countSum=graph.del[anchor].weightSum=decoded[pos+2];
			}
			for(int c=0; c<alphabet.length; c++){node.count[alphabet[c]]=node.weight[alphabet[c]]=decoded[pos+3+c];}
			node.countSum=node.weightSum=decoded[pos+1]; tail=node;
		}
		return new Family(graph, parsed-start, timed ? System.nanoTime()-parsed : 0, start-splitStart);
	}

	/** Worker times are sums and overlap each other and input; they are not additive wall phases. */
	static final class Timings {long totalNanos, readNanos, parseNanos, reconstructNanos, lineNanos, headerNanos, hashNanos, splitNanos, rawIoNanos;}
	private static final class Family {
		Family(AAGraph graph_, long parse_, long reconstruct_, long split_){graph=graph_; parse=parse_; reconstruct=reconstruct_; split=split_;}
		final AAGraph graph;
		final long parse, reconstruct, split;
	}
	private static void append(Family family, ArrayList<AAGraph> graphs, Timings timings){
		assert(family.graph!=null) : "A failed compact family must throw before ordered publication";
		graphs.add(family.graph);
		if(timings!=null){timings.parseNanos+=family.parse; timings.reconstructNanos+=family.reconstruct; timings.splitNanos+=family.split;}
	}
	private static Family get(Future<Family> future) throws IOException{
		try{return future.get();}
		catch(InterruptedException e){Thread.currentThread().interrupt(); throw new IOException("Interrupted compact HBM load", e);}
		catch(ExecutionException e){
			final Throwable cause=e.getCause();
			if(cause instanceof Error){throw (Error)cause;}
			if(cause instanceof RuntimeException){throw (RuntimeException)cause;}
			throw new IOException("Compact HBM worker failed", cause);
		}
	}

	static HbmBundleLoader.Loaded loadVerified(String path, List<String> roster,
			HbmBundleLoader.ConsensusProvider consensus, byte[][] trusted, int minCount) throws IOException{
		if(roster==null || consensus==null || trusted==null || (trusted.length!=8 && trusted.length!=16)){
			throw bad("roster, consensus and8/16 trusted semantic pins required");
		}
		final Result result=read(path);
		if(roster.size()!=result.ids.length){throw bad("external roster size differs");}
		for(int i=0; i<trusted.length; i++){
			if(trusted[i]==null || !result.provenance[i].equals(DigestSuffix.fromDigest(trusted[i]))){throw bad("semantic provenance mismatch "+i);}
		}
		for(int i=0; i<result.ids.length; i++){
			if(!result.ids[i].equals(roster.get(i)) || !Arrays.equals(result.graphs[i].pivot, consensus.consensusFor(result.ids[i]))){throw bad("external family/consensus differs");}
			if(minCount>1){
				final AAGraph graph=result.graphs[i];
				for(int p=0; p<graph.ref.length; p++){
					HbmDenseTextLoader.filterCounts(graph.ref[p], minCount, graph.pivot[p]);
					HbmDenseTextLoader.filterChain(graph.ref[p], minCount);
				}
			}
		}
		return result.loaded;
	}

	/** Graphs are privately owned; exposed consensus arrays are copies. */
	public static final class Result {
		Result(String[] ids_, AAGraph[] graphs_, HbmCompactMetadata[] metadata_, String[] provenance_, String filter_){
			ids=ids_; graphs=graphs_; metadata=metadata_; provenance=provenance_; filter=filter_;
			loaded=new HbmBundleLoader.Loaded(ids, graphs);
		}
		public int size(){return ids.length;}
		public String name(int rank){return ids[rank];}
		public byte[] consensus(int rank){return graphs[rank].pivot.clone();}
		public HbmCompactMetadata metadata(int rank){return metadata[rank];}
		public HbmBundleLoader.Loaded models(){return loaded;}
		public final String filter;
		private final String[] ids, provenance;
		private final AAGraph[] graphs;
		private final HbmCompactMetadata[] metadata;
		private final HbmBundleLoader.Loaded loaded;
	}

	/** Writer and reader use the same normalized bytes; checksum is not a model-author signature. */
	static final class Output {
		Output(String path){writer=new ByteStreamWriter(path, false, false, false); writer.start();}
		void print(String text){print(new ByteBuilder(text));}
		void print(ByteBuilder text){digest.update(text.array, 0, text.length()); writer.print(text);}
		void finish(){writer.print("#checksum\t"+DigestSuffix.fromDigest(digest.digest())+"\n");}
		boolean close(){return writer.poisonAndWait();}
		final ByteStreamWriter writer;
		final MessageDigest digest=HbmDenseTextBundle.digest();
	}

	private static final class Input {
		Input(String path, boolean blocks, Timings timings){
			file=blocks ? null : ByteFile.makeByteFile(path, false);
			raw=blocks ? new BlockInput(path, timings!=null) : null; timed=timings!=null;
		}
		byte[] next() throws IOException{
			final long start=timed ? System.nanoTime() : 0;
			final byte[] line=raw==null ? file.nextLine() : raw.next();
			if(timed){lineNanos+=System.nanoTime()-start;}
			if(line!=null){
				if(line.length==0 || (bytes+=line.length+1L)>HbmBundleFormat.CAP_TOTAL_BYTES){throw bad("empty line or text byte cap exceeded");}
				if(!starts(line, "#checksum\t")){
					final long hashStart=timed ? System.nanoTime() : 0;
					HbmDenseTextBundle.update(digest, line);
					if(timed){hashNanos+=System.nanoTime()-hashStart;}
				}
			}
			return line;
		}
		Rows rows(int count) throws IOException{
			if(raw==null){
				final ArrayList<byte[]> lines=new ArrayList<byte[]>(); long size=0;
				for(int r=0; r<count; r++){
					final byte[] line=required(next()); size+=line.length+1L;
					if(size>HbmBundleLoader.CAP_BLOCK_BYTES){throw bad("family text exceeds byte cap");}
					lines.add(line);
				}
				return new Rows(lines, null, count);
			}
			final long start=timed ? System.nanoTime() : 0;
			final byte[] block=raw.rows(count);
			if(timed){lineNanos+=System.nanoTime()-start;}
			if((bytes+=block.length)>HbmBundleFormat.CAP_TOTAL_BYTES){throw bad("text byte cap exceeded");}
			final long hashStart=timed ? System.nanoTime() : 0;
			digest.update(block);
			if(timed){hashNanos+=System.nanoTime()-hashStart;}
			return new Rows(null, block, count);
		}
		boolean close(){return file==null ? raw.close() : file.close();}
		long bytes, lineNanos, hashNanos;
		final boolean timed;
		final ByteFile file;
		final BlockInput raw;
		final LineParser1 fields=new LineParser1('\t');
		final MessageDigest digest=HbmDenseTextBundle.digest();
	}

	/** Body lines are split by the owning family worker, never by the framing reader. */
	private static final class Rows {
		Rows(ArrayList<byte[]> lines_, byte[] block_, int count_){lines=lines_; block=block_; count=count_;}
		ArrayList<byte[]> lines(){
			if(lines!=null){return lines;}
			final ArrayList<byte[]> result=new ArrayList<byte[]>(); int start=0;
			for(int i=0; i<block.length; i++){
				if(block[i]=='\n'){
					if(i==start){throw bad("empty body line");}
					result.add(Arrays.copyOfRange(block, start, i)); start=i+1;
				}
			}
			if(start!=block.length || result.size()!=count){throw bad("body row count or termination mismatch");}
			return result;
		}
		final ArrayList<byte[]> lines;
		final byte[] block;
		final int count;
	}

	/** Native decompression plus bounded #rows framing. CRLF bodies normalize before hashing. */
	private static final class BlockInput {
		BlockInput(String path_, boolean timed_){path=path_; timed=timed_; input=ReadWrite.getInputStream(path, false, false);}
		boolean fill() throws IOException{
			if(pos<limit){return true;}
			final long start=timed ? System.nanoTime() : 0;
			do{limit=input.read(buffer);}while(limit==0);
			if(timed){ioNanos+=System.nanoTime()-start;}
			pos=0; lineIndex=0; newlines.clear();
			if(limit>0){
				// The Java8-safe Vector facade retains the scalar fallback automatically.
				simd.Vector.findSymbols(buffer, 0, limit, (byte)'\n', newlines);
				hasCR=simd.Vector.countSymbols(buffer, 0, limit, (byte)'\r')>0;
			}
			return limit>0;
		}
		byte[] next() throws IOException{
			final ByteBuilder line=new ByteBuilder();
			while(fill()){
				final int start=pos; while(pos<limit && buffer[pos]!='\n'){pos++;}
				append(line, start, pos-start);
				if(pos<limit){
					pos++; lineIndex++; if(line.length()>0 && line.array[line.length()-1]=='\r'){line.setLength(line.length()-1);}
					return line.toBytes();
				}
			}
			return line.length()==0 ? null : line.toBytes();
		}
		byte[] rows(int count) throws IOException{
			assert(count>0) : "Compact #rows is validated against a nonempty consensus before framing";
			final ByteBuilder block=new ByteBuilder(); int rows=0; boolean cr=false;
			while(rows<count){
				if(!fill()){throw new IOException("Truncated compact HBM body");}
				final int start=pos;
				final int take=Math.min(count-rows, newlines.size()-lineIndex);
				rows+=take; lineIndex+=take; cr|=hasCR;
				pos=rows==count ? newlines.get(lineIndex-1)+1 : limit;
				append(block, start, pos-start);
			}
			if(cr){
				int length=0;
				for(int i=0; i<block.length(); i++){
					if(block.array[i]!='\r' || i+1==block.length() || block.array[i+1]!='\n'){block.array[length++]=block.array[i];}
				}
				block.setLength(length);
			}
			return block.toBytes();
		}
		void append(ByteBuilder out, int start, int length){
			if(out.length()+(long)length>HbmBundleLoader.CAP_BLOCK_BYTES){throw bad("family text exceeds byte cap");}
			out.append(buffer, start, length);
		}
		boolean close(){return ReadWrite.finishReading(input, path, false);}
		final String path;
		final boolean timed;
		final java.io.InputStream input;
		final byte[] buffer=new byte[65536];
		final structures.IntList newlines=new structures.IntList();
		int pos, limit, lineIndex;
		boolean hasCR;
		long ioNanos;
	}

	private static int[] alphabet(String text){
		if(text.isEmpty() || text.length()>22){throw bad("invalid alphabet length");}
		final int[] out=new int[text.length()]; final boolean[] seen=new boolean[22];
		for(int i=0; i<out.length; i++){
			int code=-1;
			for(int j=0; j<22; j++){if(text.charAt(i)==(j==Blosum62.X_CODE ? 'X' : (char)AminoAcid.numberToAcid[j])){code=j; break;}}
			if(code<0 || seen[code]){throw bad("unknown or duplicate alphabet symbol");}
			seen[code]=true; out[i]=code;
		}
		return out;
	}
	private static String field(LineParser1 row){if(row.terms()!=2){throw bad("expected one header value");} return row.parseString(1);}
	private static boolean starts(byte[] line, String prefix){
		if(line.length<prefix.length()){return false;}
		for(int i=0; i<prefix.length(); i++){if(line[i]!=(byte)prefix.charAt(i)){return false;}}
		return true;
	}
	private static byte[] required(byte[] line) throws IOException{if(line==null){throw new IOException("Truncated compact HBM");} return line;}
	private static void require(byte[] line, String text) throws IOException{
		if(!Arrays.equals(required(line), text.getBytes(java.nio.charset.StandardCharsets.US_ASCII))){throw bad("expected "+text);}
	}
	private static IllegalArgumentException bad(String reason){return new IllegalArgumentException("Compact HBM: "+reason);}
	private HbmCompactTextReader(){}
}
