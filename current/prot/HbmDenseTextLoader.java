package prot;

import java.io.IOException;
import java.security.MessageDigest;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import fileIO.ByteFile;
import parse.LineParser1;

/**
	 * Reader for dense, explicit, candidate-only HBM text.
 * Numeric fields are parsed directly from bytes; independent families may be
 * reconstructed concurrently. At most twice the worker count is queued.
 *
	 * The v1 dump entrypoint is experimental and has no integrity contract. The
	 * v2 entrypoint verifies semantic provenance, family checksums and the root
	 * checksum before publishing any models; see HbmDenseTextBundle.
 *
 * @author Collei
 */
public final class HbmDenseTextLoader {

	private HbmDenseTextLoader(){}

	/** Loads the strict mincount=1 subset, checking roster and consensus binding. */
	static HbmBundleLoader.Loaded load(final String path, final List<String> roster,
			final HbmBundleLoader.ConsensusProvider consensus, final String sourceSha80,
			final int threads) throws IOException{
		DigestSuffix.requireSuffix(sourceSha80, "candidate source pin");
		return loadInternal(path, roster, consensus, sourceSha80, threads, null, 1);
	}

	/** Verified v2 entrypoint. The caller supplies the same trusted pins as MQHB. */
	static HbmBundleLoader.Loaded loadVerified(final String path, final List<String> roster,
			final HbmBundleLoader.ConsensusProvider consensus, final byte[][] trusted,
			final int threads) throws IOException{
		return loadVerified(path, roster, consensus, trusted, threads, 1);
	}

	/** Applies the same optional in-memory count filter as the native reader. */
	static HbmBundleLoader.Loaded loadVerified(final String path, final List<String> roster,
			final HbmBundleLoader.ConsensusProvider consensus, final byte[][] trusted,
			final int threads, final int minCount) throws IOException{
		if(minCount<1){throw bad("minCount must be at least1");}
		if(trusted==null || (trusted.length!=8 && trusted.length!=16)){throw bad("8 or16 trusted semantic pins required");}
		for(final byte[] pin : trusted){
			if(pin==null || (pin.length!=10 && pin.length!=32)){throw bad("invalid trusted semantic pin width");}
		}
		return loadInternal(path, roster, consensus, null, threads, trusted, minCount);
	}

	private static HbmBundleLoader.Loaded loadInternal(final String path, final List<String> roster,
			final HbmBundleLoader.ConsensusProvider consensus, final String sourceSha80,
			final int threads, final byte[][] trusted, final int minCount) throws IOException{
		if(threads<1 || threads>256){throw new IllegalArgumentException("threads must be in [1,256]");}
		if(roster==null || roster.isEmpty() || consensus==null){throw new IllegalArgumentException("roster and consensus required");}
		final HashSet<String> seen=new HashSet<String>();
		for(final String id : roster){
			if(id==null || id.isEmpty() || !seen.add(id) || HbmBundleLoader.utf8StrictEncode(id).length>HbmBundleFormat.CAP_REPID){
				throw new IllegalArgumentException("invalid or duplicate roster ID");
			}
		}
		final String[] ids=roster.toArray(new String[0]);
		final AAGraph[] graphs=new AAGraph[ids.length];
		final ByteFile bf=ByteFile.makeByteFile(path, false);
		final ExecutorService pool=(threads==1 ? null : Executors.newFixedThreadPool(threads));
		final ArrayDeque<Future<AAGraph>> pending=new ArrayDeque<Future<AAGraph>>();
		boolean completed=false;
		final MessageDigest root=(trusted==null ? null : HbmDenseTextBundle.digest());
		try{
			header(bf, root, "#format", root==null ? "hbm_text_v1" : HbmDenseTextBundle.FORMAT);
			if(root!=null){
				header(bf, root, "#graph_contract", HbmDenseTextBundle.CONTRACT);
				for(int i=0; i<16; i++){
					final String observed=header(bf, root, "#provenance_"+i, null);
					DigestSuffix.requireSuffix(observed, "provenance "+i);
					if(i<trusted.length && !observed.equals(DigestSuffix.fromDigest(trusted[i]))){throw bad("semantic provenance mismatch at "+i);}
					if((i==1 || i==2) && observed.equals("00000000000000000000")){throw bad("zero loader provenance");}
				}
			}
			final String source=header(bf, root, "#candidate_sha80", sourceSha80);
			DigestSuffix.requireSuffix(source, "candidate source pin");
			header(bf, root, "#original_sha80", source);
			header(bf, root, "#sum_source", "candidate");
			header(bf, root, "#n_families", Integer.toString(ids.length));
			header(bf, root, "#encoding", "dense");
			header(bf, root, "#coordinate_mode", "explicit");
			header(bf, root, "#min_count_emitted", "1");
			header(bf, root, "#dense_row_threshold_nonzero_symbols", "11");
			header(bf, root, "#start_end_counts", "unavailable");
			// Descriptive rows emitted by HbmTextDump; reject missing/reordered headers.
			for(final String key : new String[]{"#row_f", "#row_r", "#row_i", "#row_d", "#row_di", "#symbol_count"}){
				final LineParser1 lp=new LineParser1('\t');
				final byte[] line=required(bf);
				lp.set(line);
				if(lp.terms()!=2 || !lp.termEquals(key, 0)){throw bad("missing descriptive header "+key);}
				HbmDenseTextBundle.update(root, line);
			}
			long bytes=0, nodes=0;
			int joined=0;
			for(int rank=0; rank<ids.length; rank++){
				final ArrayList<byte[]> lines=new ArrayList<byte[]>();
				final byte[] first=required(bf);
				final LineParser1 lp=new LineParser1('\t');
				lp.set(first);
				if(lp.terms()!=6 || !lp.termEquals('f', 0) || integer(lp, 1)!=rank){throw bad("family rank/header at "+rank);}
				final int length=integer(lp, 3);
				if(length<1 || length>HbmBundleFormat.CAP_L){throw bad("family length "+length);}
				nodes+=2L*length;
				if(nodes>HbmBundleFormat.CAP_TOTAL_NODES){throw bad("total node cap exceeded");}
				lines.add(first);
				long familyBytes=first.length+1L;
				byte[] line=required(bf);
				for(; !(line.length>0 && line[0]=='e'); line=required(bf)){
					familyBytes+=line.length+1L;
					if(familyBytes>HbmBundleLoader.CAP_BLOCK_BYTES){throw bad("family text exceeds byte cap");}
					if(line.length>0 && line[0]=='i' || line.length>1 && line[0]=='d' && line[1]=='i'){
						if(++nodes>HbmBundleFormat.CAP_TOTAL_NODES){throw bad("total node cap exceeded");}
					}
					lines.add(line);
				}
				final String familyPin;
				if(root==null){
					if(line.length!=1){throw bad("invalid v1 family terminator");}
					familyPin=null;
				}else{
					lp.set(line);
					if(lp.terms()!=2 || !lp.termEquals('e', 0)){throw bad("invalid v2 family checksum row");}
					familyPin=DigestSuffix.requireSuffix(lp.parseString(1), "family checksum");
					HbmDenseTextBundle.update(root, line);
				}
				bytes+=familyBytes+2;
				if(bytes>HbmBundleFormat.CAP_TOTAL_BYTES){throw bad("total text exceeds byte cap");}
				final String id=ids[rank];
				final byte[] provided=consensus.consensusFor(id);
				if(provided==null || provided.length!=length){throw bad("missing/wrong consensus for "+id);}
				final byte[] pivot=provided.clone();
				if(pool==null){graphs[rank]=family(lines, id, pivot, familyPin, minCount);}
				else{
					pending.addLast(pool.submit(() -> family(lines, id, pivot, familyPin, minCount)));
					if(pending.size()>=2*threads){graphs[joined++]=get(pending.removeFirst());}
				}
			}
			final byte[] end=required(bf);
			if(root==null){
				if(end.length!=1 || end[0]!='z'){throw bad("missing final z");}
			}else{
				final LineParser1 lp=new LineParser1('\t');
				lp.set(end);
				if(lp.terms()!=2 || !lp.termEquals('z', 0) || !lp.termEquals(DigestSuffix.fromDigest(root.digest()), 1)){
					throw bad("root checksum mismatch");
				}
			}
			if(bf.nextLine()!=null){throw bad("trailing data after final z");}
			while(!pending.isEmpty()){graphs[joined++]=get(pending.removeFirst());}
			completed=true;
			return new HbmBundleLoader.Loaded(ids, graphs);
		}finally{
			if(pool!=null){
				pool.shutdownNow();
				boolean interrupted=false;
				for(;;){
					try{if(pool.awaitTermination(1, TimeUnit.DAYS)){break;}}
					catch(InterruptedException e){interrupted=true;}
				}
				if(interrupted){Thread.currentThread().interrupt();}
			}
			final boolean error=bf.close();
			if(error && completed){throw new IOException("Failed closing dense text input: "+path);}
		}
	}

	/** Rebuilds one family's two insertion namespaces and all retained counts. */
	private static AAGraph family(final ArrayList<byte[]> lines, final String id, final byte[] pivot,
			final String pin, final int minCount){
		if(pin!=null){
			final MessageDigest md=HbmDenseTextBundle.digest();
			for(final byte[] line : lines){HbmDenseTextBundle.update(md, line);}
			if(!pin.equals(DigestSuffix.fromDigest(md.digest()))){throw bad("family checksum mismatch: "+id);}
		}
		final LineParser1 lp=new LineParser1('\t');
		lp.set(lines.get(0));
		final String encoded=HbmBundleLoader.utf8StrictDecode(lp.parseByteArray(2));
		if(!id.equals(IdentifierCodec.decode(encoded)) || !IdentifierCodec.encode(id).equals(encoded)){throw bad("noncanonical/wrong family ID: "+id);}
		final byte[] raw=lp.parseByteArray(5);
		final byte[] observed=Blosum62.encode(raw, id);
		if(!Arrays.equals(pivot, observed) || !lp.termEquals(DigestSuffix.bytes(pivot), 4)){throw bad("consensus mismatch: "+id);}
		for(final byte b : pivot){if(!((b>=0 && b<20) || b==Blosum62.X_CODE)){throw bad("invalid consensus residue");}}
		final AAGraph graph=new AAGraph(pivot, 0);
		HbmBundleLoader.checkGraphKnobs(graph);
		int anchor=-1, refChain=0, delChain=0, phase=0;
		AAGraphNode refTail=null, delTail=null;
		for(int row=1; row<lines.size(); row++){
			if(Thread.currentThread().isInterrupted()){throw bad("family parsing interrupted");}
			lp.set(lines.get(row));
			if(lp.terms()<3){throw bad("short family row in "+id);}
			final int pos=integer(lp, 1);
			if(lp.termEquals('r', 0)){
				if(pos!=anchor+1 || pos>=pivot.length){throw bad("REF order/range in "+id);}
				anchor=pos; refChain=delChain=phase=0;
				refTail=graph.ref[pos]; delTail=graph.del[pos];
				counts(lp, refTail, 2);
				if(refTail.count[pivot[pos]]<1){throw bad("missing pivot count in "+id);}
			}else{
				if(anchor<0 || pos!=anchor){throw bad("row precedes its REF anchor in "+id);}
				if(lp.termEquals('d', 0)){
					if(lp.terms()!=3 || phase!=0){throw bad("duplicate or misplaced DEL in "+id);}
					phase=1;
					final int sum=integer(lp, 2);
					if(sum<1){throw bad("zero DEL must be omitted");}
					delTail.countSum=delTail.weightSum=sum;
				}else if(lp.termEquals('i', 0) || lp.termEquals("di", 0)){
					final boolean del=lp.termEquals("di", 0);
					final int index=integer(lp, 2), expected=(del ? delChain : refChain);
					if(index!=expected || index>=HbmBundleFormat.CAP_CHAIN || (!del && phase!=0)){throw bad("insertion chain order/cap in "+id);}
					final AAGraphNode node=new AAGraphNode(Blosum62.X_CODE, AAGraphNode.INS, anchor+1);
					counts(lp, node, 3);
					if(del){phase=2; delChain++; delTail.insEdge=node; delTail=node;}
					else{refChain++; refTail.insEdge=node; refTail=node;}
				}else{throw bad("unknown row in "+id);}
			}
		}
		if(anchor!=pivot.length-1){throw bad("missing final REF rows in "+id);}
		if(minCount>1){
			// Validate the complete stored graph before filtering, just as the native
			// reader validates every stored count even beyond a truncated INS chain.
			for(int i=0; i<pivot.length; i++){
				filterCounts(graph.ref[i], minCount, pivot[i]);
				filterChain(graph.ref[i], minCount);
				filterChain(graph.del[i], minCount);
			}
		}
		return graph;
	}

	/** Mirrors HbmBundleLoader.overwrite: protect REF pivot, retain DEL occupancy. */
	static void filterCounts(final AAGraphNode node, final int minCount, final int protectedResidue){
		assert(node.type!=AAGraphNode.DEL) : "DEL has no histogram; native minCount only filters REF/INS slots";
		int sum=0;
		for(int i=0; i<node.count.length; i++){
			final int retained=(node.count[i]>=minCount || i==protectedResidue) ? node.count[i] : 0;
			node.count[i]=node.weight[i]=retained;
			sum=Math.addExact(sum, retained);
		}
		node.countSum=node.weightSum=sum;
	}

	/** An empty filtered INS node truncates its whole suffix, per native readChain. */
	static void filterChain(final AAGraphNode parent, final int minCount){
		AAGraphNode previous=parent;
		for(AAGraphNode node=parent.insEdge; node!=null; node=node.insEdge){
			filterCounts(node, minCount, -1);
			if(node.countSum==0){previous.insEdge=null; break;}
			previous=node;
		}
	}

	/** Dense histograms have native NAA slots; native builder requires weight=count. */
	private static void counts(final LineParser1 lp, final AAGraphNode node, final int mode){
		if(lp.terms()!=mode+2+HbmBundleFormat.NAA || !lp.termEquals('d', mode)){throw bad("expected dense 22-slot histogram");}
		final int sum=integer(lp, mode+1);
		long actual=0;
		for(int i=0; i<HbmBundleFormat.NAA; i++){
			final int count=integer(lp, mode+2+i);
			node.count[i]=node.weight[i]=count;
			actual+=count;
		}
		if(sum<1 || actual!=sum){throw bad("histogram sum mismatch: "+actual+" != "+sum);}
		node.countSum=node.weightSum=sum;
	}

	/** Canonical, overflow-checked A48 int parsing without per-field allocations. */
	static int integer(final LineParser1 lp, final int term){
		return integer(lp, term, 0);
	}

	/** Compact DEL tokens prefix the same canonical integer with one '-'. */
	static int integer(final LineParser1 lp, final int term, final int skip){
		final int length=lp.length(term)-skip, start=lp.a()+skip, end=lp.b();
		final byte[] line=lp.line();
		if(length<1 || length>6 || (length>1 && line[start]=='0')){throw bad("noncanonical A48 integer at field "+term);}
		long value=0;
		for(int i=start; i<end; i++){
			final int digit=line[i]-48;
			if(digit<0 || digit>63){throw bad("invalid A48 digit at field "+term);}
			value=(value<<6)|digit;
		}
		if(value>Integer.MAX_VALUE){throw bad("A48 integer overflow at field "+term);}
		return (int)value;
	}

	private static String header(final ByteFile bf, final MessageDigest root, final String key, final String value) throws IOException{
		final LineParser1 lp=new LineParser1('\t');
		final byte[] line=required(bf);
		lp.set(line);
		if(lp.terms()!=2 || !lp.termEquals(key, 0) || (value!=null && !lp.termEquals(value, 1))){throw bad("expected header "+key+"="+value);}
		HbmDenseTextBundle.update(root, line);
		return value==null ? lp.parseString(1) : value;
	}
	private static byte[] required(final ByteFile bf) throws IOException{
		final byte[] line=bf.nextLine();
		if(line==null){throw new IOException("Unexpected end of dense HBM text");}
		return line;
	}
	private static AAGraph get(final Future<AAGraph> future) throws IOException{
		try{return future.get();}
		catch(InterruptedException e){Thread.currentThread().interrupt(); throw new IOException("Interrupted dense HBM load", e);}
		catch(ExecutionException e){
			final Throwable cause=e.getCause();
			if(cause instanceof Error){throw (Error)cause;}
			if(cause instanceof RuntimeException){throw (RuntimeException)cause;}
			throw new IOException("Dense HBM worker failed", cause);
		}
	}
	private static IllegalArgumentException bad(final String message){return new IllegalArgumentException("HBM dense text: "+message);}
}
