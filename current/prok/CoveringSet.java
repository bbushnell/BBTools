package prok;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import dna.AminoAcid;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import map.LongHashSet;
import map.LongIntMap;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prot.ReducedAlphabet;
import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.FASTQ;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import stream.Writer;
import stream.WriterFactory;
import structures.LongList;
import structures.ListNum;

/**
 * Greedy set-cover kmer selection for tRNA (or any short-gene) covering sets.
 * Replaces the bash kmercountexact+bbduk descent loop with a single in-memory
 * pass: count all kmers across a pool of sequences, then iteratively select the
 * top-N by count, evict covered sequences, recount, and repeat until a coverage
 * target or kmer budget is reached.
 *
 * Two-tier ranking (Brian's design): each round takes the top 2*step kmers by
 * CURRENT count on the remaining pool, re-sorts them by ORIGINAL prevalence
 * (recorded before any eviction), and keeps the top step. This concentrates on
 * kmers shared between rare and common tRNAs rather than random sequence that
 * happens to be common among the rare stragglers.
 *
 * @author Neptune, Brian Bushnell
 */
public class CoveringSet {

	public static void main(String[] args){
		if(args.length==1 && args[0].equalsIgnoreCase("selftest")){
			selftest();
			return;
		}
		Timer t=new Timer();
		CoveringSet x=new CoveringSet(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	/** Small paper-checkable regression for the protein/family path. */
	private static void selftest(){
		try{
			final java.nio.file.Path dir=java.nio.file.Files.createTempDirectory("covering-set-selftest.");
			final java.nio.file.Path families=java.nio.file.Files.createDirectory(dir.resolve("families"));
			final java.nio.file.Path a=families.resolve("A.faa"), b=families.resolve("B.faa");
			final java.nio.file.Path exclude=dir.resolve("exclude.tsv");
			java.nio.file.Files.write(a, ">a1\nACDEFG\n>a2\nACDEFW\n>a3\nWWWWAC\n".getBytes(java.nio.charset.StandardCharsets.US_ASCII));
			java.nio.file.Files.write(b, ">b1\nACDEFG\n>b2\nHHHHAC\n>b3\nGGGGGG\n".getBytes(java.nio.charset.StandardCharsets.US_ASCII));
			java.nio.file.Files.write(exclude, "#id\taction\na2\tholdout\n".getBytes(java.nio.charset.StandardCharsets.US_ASCII));
			final java.nio.file.Path aaOut=dir.resolve("aa20.sets.tsv"), aaSummary=dir.resolve("aa20.summary.tsv");
			new CoveringSet(new String[]{"families="+families, "out="+aaOut, "summary="+aaSummary,
				"alphabet=amino", "k=2", "step=1", "target=1.0", "minhits=2", "t=1"}).process(new Timer());
			final java.nio.file.Path c8Out=dir.resolve("c8.sets.tsv"), c8Summary=dir.resolve("c8.summary.tsv");
			new CoveringSet(new String[]{"families="+families, "out="+c8Out, "summary="+c8Summary,
				"alphabet=ACDFGHIW", "key=AST/C/DEKNQR/FY/GP/H/ILMV/W", "exclude="+exclude,
				"maxfamilies=1", "k=2", "step=1", "target=1.0", "minhits=2", "t=1"}).process(new Timer());
			final java.nio.file.Path c8T32Out=dir.resolve("c8.t32.sets.tsv"), c8T32Summary=dir.resolve("c8.t32.summary.tsv");
			new CoveringSet(new String[]{"families="+families, "out="+c8T32Out, "summary="+c8T32Summary,
				"alphabet=ACDFGHIW", "key=AST/C/DEKNQR/FY/GP/H/ILMV/W", "exclude="+exclude,
				"maxfamilies=1", "k=2", "step=1", "target=1.0", "minhits=2", "t=32"}).process(new Timer());
			final String aa=new String(java.nio.file.Files.readAllBytes(aaOut), java.nio.charset.StandardCharsets.US_ASCII);
			final String c8=new String(java.nio.file.Files.readAllBytes(c8Out), java.nio.charset.StandardCharsets.US_ASCII);
			final String sum=new String(java.nio.file.Files.readAllBytes(c8Summary), java.nio.charset.StandardCharsets.US_ASCII);
			final String c8T32=new String(java.nio.file.Files.readAllBytes(c8T32Out), java.nio.charset.StandardCharsets.US_ASCII);
			final String sumT32=new String(java.nio.file.Files.readAllBytes(c8T32Summary), java.nio.charset.StandardCharsets.US_ASCII);
			if(!aa.contains("#columns\tfamily_id\tkmer") || !c8.contains("#columns\tfamily_id\tkmer")){throw new RuntimeException("missing covering-set output schema");}
			if(!sum.contains("A\t3\t1\t") || !sum.contains("B\t3\t0\t")){throw new RuntimeException("exclude/member accounting failed");}
			if(!c8.contains("#cross_family_candidates_removed\t")){throw new RuntimeException("maxfamilies accounting missing");}
			if(!c8.equals(c8T32) || !sum.equals(sumT32)){throw new RuntimeException("t=1 and t=32 outputs differ");}
			System.out.println("PASS CoveringSet selftest: aa20+c8, six-sequence pool, exclude, minhits=2, maxfamilies=1, deterministic t=1/32");
		} catch(Exception e){throw new RuntimeException("CoveringSet selftest failed", e);}
	}

	public CoveringSet(String[] args){
		{
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		FASTQ.TEST_INTERLEAVED=FASTQ.FORCE_INTERLEAVED=false;

		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("k")){
				k=Integer.parseInt(b);
			}else if(a.equals("kdesign")){
				kDesign=Integer.parseInt(b);
			}else if(a.equals("step")){
				step=(int)Parse.parseKMG(b);
			}else if(a.equals("stepfraction") || a.equals("stepfrac") || a.equals("adaptivefraction")){
				stepFraction=Double.parseDouble(b);
			}else if(a.equals("minstepmult")){
				minStepMult=Double.parseDouble(b);
			}else if(a.equals("maxstepmult")){
				maxStepMult=Double.parseDouble(b);
			}else if(a.equals("step2boost")){
				step2Boost=Double.parseDouble(b);
			}else if(a.equals("maxkmers") || a.equals("budget")){
				maxKmers=(int)Parse.parseKMG(b);
			}else if(a.equals("mincov") || a.equals("mincoverage") || a.equals("target")){
				minCovFraction=Float.parseFloat(b);
			}else if(a.equals("extra") || a.equals("extra1")){
				extra=b;
			}else if(a.equals("copies")){
				copies=Integer.parseInt(b);
			}else if(a.equals("rcomp")){
				rcomp=Parse.parseBoolean(b);
			}else if(a.equals("alphabet")){
				alphabetSpec=b;
			}else if(a.equals("key")){
				keySpec=b;
			}else if(a.equals("minhits")){
				minHits=Integer.parseInt(b);
			}else if(a.equals("families") || a.equals("familydir")){
				families=b;
			}else if(a.equals("exclude")){
				exclude=b;
			}else if(a.equals("summary")){
				summary=b;
			}else if(a.equals("maxfamilies")){
				maxFamilies=Integer.parseInt(b);
			}else if(a.equals("bbtools_commit") || a.equals("commit")){
				bbtoolsCommit=b;
			}else if(a.equals("partitions") || a.equals("numpartitions")){
				numPartitions=Integer.parseInt(b);
			}else if(a.equals("bufsize") || a.equals("partitionbuffer")){
				bufferSize=Integer.parseInt(b);
			}else if(parser.parse(arg, a, b)){
				//handled
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}

		in=parser.in1;
		out=parser.out1;
		if(families==null && in==null){throw new RuntimeException("Error - an input file is required.");}
		if(out==null){throw new RuntimeException("Error - an output file is required.");}
		if(families!=null && summary==null){throw new RuntimeException("families= requires summary=");}
		if(families!=null && out.equals(summary)){throw new IllegalArgumentException("out= and summary= must be different files");}
		if(keySpec!=null && (alphabetSpec==null || alphabetSpec.equalsIgnoreCase("nt"))){
			throw new IllegalArgumentException("key= requires an amino alphabet");
		}
		if(minHits<1){throw new IllegalArgumentException("minhits must be >=1: "+minHits);}
		if(maxFamilies<0){throw new IllegalArgumentException("maxfamilies must be >=0: "+maxFamilies);}
		if(alphabetSpec==null || alphabetSpec.equalsIgnoreCase("nt")){
			proteinMode=false;
		}else{
			reducedAlphabet=ReducedAlphabet.parse(alphabetSpec, keySpec);
			proteinMode=true;
			if(rcomp){throw new IllegalArgumentException("rcomp is valid only for alphabet=nt");}
		}
		if(families!=null && !proteinMode){throw new IllegalArgumentException("families= requires an amino alphabet");}

		if(kDesign<0){kDesign=k;}
		if(k<1 || k>31){throw new IllegalArgumentException("k must be 1..31: "+k);}
		if(kDesign<k || kDesign>31){throw new IllegalArgumentException("kdesign must be in [k,31]: "+kDesign);}
		if(proteinMode && (long)kDesign*reducedAlphabet.bits()>62){
			throw new IllegalArgumentException("Packed k-mer exceeds 62 bits: k="+kDesign+
				" bits="+reducedAlphabet.bits()+" alphabet="+reducedAlphabet.symbols());
		}
		if(step<1){throw new IllegalArgumentException("step must be positive: "+step);}
		if(stepFraction<0 || stepFraction>1){
			throw new IllegalArgumentException("stepfraction must be in [0,1]: "+stepFraction);
		}
		if(minStepMult<=0 || maxStepMult<minStepMult || step2Boost<=0){
			throw new IllegalArgumentException("Require 0<minstepmult<=maxstepmult and step2boost>0: "+
					minStepMult+", "+maxStepMult+", "+step2Boost);
		}

		overwrite=parser.overwrite;
		ffin=(in==null ? null : FileFormat.testInput(in, FileFormat.FASTA, null, true, true));
		ffout=(families!=null ? null : FileFormat.testOutput(out, FileFormat.FASTA, null, true, overwrite, false, false));
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	void process(Timer t){
		if(families!=null){processFamilies(t); return;}
		if(proteinMode){processSingleProtein(t); return;}
		processNucleotide(t);
	}

	/** Original nucleotide implementation.  This branch is intentionally kept
	 * separate so alphabet extensions cannot perturb the accepted nt bytes. */
	private void processNucleotide(Timer t){
		ArrayList<byte[]> pool=loadPool();
		final int totalSeqs=pool.size();
		//numPartitions unset (<=0) means "scale with available threads" -- resolved
		//here rather than in field init since Shared.threads() reflects t= parsing
		//that happens after this object is constructed. Default is ~2x the thread
		//count (extra shards make it less likely two threads want to flush the same
		//partition at once, reducing lock contention) and forced ODD, since way()'s
		//kmer%numPartitions would otherwise alias against the structured low-bit
		//patterns of 2-bit-packed kmers -- an even modulus can route every kmer
		//ending in the same base to the same parity of partitions (Brian, 2026-08-28).
		if(numPartitions<=0){numPartitions=Tools.max(15, 2*Shared.threads())|1;}
		outstream.println("Pool: "+totalSeqs+" sequences, k="+k+
				(kDesign!=k ? " (design k="+kDesign+")" : "")+
				", step="+step+(stepFraction>0 ? ", stepfraction="+stepFraction+
				", stepbounds="+minStepMult+"-"+maxStepMult+", step2boost="+step2Boost : "")+
				", rcomp="+rcomp+", partitions="+numPartitions+", bufsize="+bufferSize);

		long totalBases=0;
		for(byte[] seq : pool){totalBases+=seq.length;}

		final long designMask=~((-1L)<<(2*kDesign));
		//totalBases is a safe (if loose) upper bound on distinct kDesign-mers -- used
		//only to pre-size the partition shards, never a correctness constraint.
		final PartitionedCounter originalCounts=countKmersMT(pool, null, kDesign, designMask, totalBases);
		outstream.println("Distinct "+kDesign+"-mers in pool: "+originalCounts.size());

		LongList selectedKmers=new LongList();
		LongHashSet selectedSet=new LongHashSet(maxKmers>0 ? maxKmers : 16384);
		boolean[] alive=new boolean[pool.size()];
		Arrays.fill(alive, true);
		int aliveCount=totalSeqs;
		int round=0;
		int currentStep=step;

		while(aliveCount>0){
			if(maxKmers>0 && selectedKmers.size>=maxKmers){break;}
			float covFrac=1f-(aliveCount/(float)totalSeqs);
			if(covFrac>=minCovFraction){break;}

			//Retain the historical near-target reduction when adaptive sizing is off.
			if(stepFraction<=0){
				final float gap=aliveCount/(float)totalSeqs;
				final float targetGap=1f-minCovFraction;
				if(gap<=2*targetGap && currentStep>Math.max(1, step/2)){
					currentStep=Math.max(1, step/2);
					outstream.println("  Step reduced to "+currentStep+" (approaching target, gap="+
							String.format("%.4f", gap)+", targetGap="+String.format("%.4f", targetGap)+")");
				}
			}
			final int roundStep=maxKmers>0 ? Tools.min(currentStep, maxKmers-selectedKmers.size) : currentStep;
			if(roundStep<=0){break;}

			//originalCounts.size() upper-bounds any later round's distinct count
			//(coverage only shrinks the alive set), so it's a safe, already-known
			//sizing hint -- no extra pool scan needed each round.
			PartitionedCounter currentCounts=countKmersMT(pool, alive, kDesign, designMask, originalCounts.size());
			if(currentCounts.isEmpty()){break;}

			final int candidateLimit=(int)Tools.min(Integer.MAX_VALUE, 2L*roundStep);
			long[] candidates=topNByCount(currentCounts, candidateLimit);
			long[] ranked=rankByOriginal(candidates, originalCounts, roundStep);

			int added=0;
			for(long kmer : ranked){
				if(selectedSet.contains(kmer)){continue;}
				selectedSet.add(kmer);
				selectedKmers.add(kmer);
				added++;
			}
			if(added==0){break;}

			int evicted=evictCoveredMT(pool, alive, selectedSet, kDesign, designMask);
			aliveCount-=evicted;
			round++;

			float cov=1f-(aliveCount/(float)totalSeqs);
			if(verbose || round%10==0){
				outstream.println("  Round "+round+": +"+added+" kmers (total "+
						selectedKmers.size+"), evicted "+evicted+
						", alive="+aliveCount+", coverage="+String.format("%.4f", cov));
			}
			if(stepFraction>0 && aliveCount>0){
				final int nextStep=adaptiveStep(added, evicted, aliveCount, round+1);
				if(verbose || nextStep!=currentStep){
					outstream.println("  Next step: "+nextStep+" (previous added="+added+
							", recovered="+evicted+", remaining="+aliveCount+")");
				}
				currentStep=nextStep;
			}
		}

		float finalCov=1f-(aliveCount/(float)totalSeqs);
		outstream.println("Selected "+selectedKmers.size+" "+kDesign+"-mers in "+round+
				" rounds, coverage="+String.format("%.6f", finalCov)+
				" ("+aliveCount+" uncovered of "+totalSeqs+")");

		ArrayList<Read> output;
		if(kDesign==k){
			output=kmersToReads(selectedKmers, k);
		}else{
			output=designToUseKmers(selectedKmers, kDesign, k);
		}
		outstream.println("Output: "+output.size()+" "+k+"-mer sequences");

		if(ffout!=null){
			//Writer, not ConcurrentReadOutputStream -- CRIS/CROS are being phased
			//out in favor of Streamer/Writer (Brian, 2026-08-26), same reasoning
			//as the loadFasta switch to Streamer above.
			Writer writer=WriterFactory.makeWriter(ffout);
			writer.start();
			writer.add(output, 0);
			errorState|=writer.poisonAndWait();
		}

		if(execPool!=null){execPool.shutdown();}

		t.stop();
		outstream.println("Time:   "+t);
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state.");
		}
	}

	/** Single-pool arbitrary-alphabet mode.  The historical nt mode above remains
	 * the byte-preserving implementation; this path is intentionally simpler and
	 * uses the same deterministic greedy ordering as the family mode. */
	private void processSingleProtein(Timer t){
		final ProteinLoad loaded=loadProteinMembers(in, Collections.<String>emptySet(), "single");
		final FamilyPool pool=new FamilyPool("single", loaded.members, in, hashFile(new File(in)), loaded.totalMembers);
		final FamilyResult result=selectFamily(pool, null);
		final PrintWriter pw=openTextOutput(out);
		for(OutputKmer row : result.output){pw.println(">kmer_"+row.rank); pw.println(row.text);}
		pw.close();
		outstream.println("Selected "+result.output.size()+" protein kmers, coverage="+result.coverage);
		t.stop(); outstream.println("Time:   "+t);
	}

	/** Processes all family pools in one JVM.  Pools are loaded before workers
	 * start so malformed input fails loudly and global family-frequency limits are
	 * computed over exactly the held-out input set. */
	private void processFamilies(Timer t){
		final Set<String> excluded=loadExclude(exclude);
		final ArrayList<FamilyInput> inputs=loadFamilyInputs(families);
		final ArrayList<FamilyPool> pools=new ArrayList<FamilyPool>(inputs.size());
		final LongIntMap familyFrequency=new LongIntMap(1024);
		for(FamilyInput input : inputs){
			final ProteinLoad loaded=loadProteinMembers(input.path, excluded, input.family);
			final ArrayList<ProteinMember> members=loaded.members;
			if(members.isEmpty()){throw new RuntimeException("Family has no non-excluded members: "+input.family);}
			final FamilyPool pool=new FamilyPool(input.family, members, input.path, input.sha256, loaded.totalMembers);
			pools.add(pool);
			final LongHashSet seen=new LongHashSet(1024);
			for(ProteinMember member : members){for(long kmer : kmers(member.bases, kDesign)){seen.add(kmer);}}
			if(seen.isEmpty()){throw new RuntimeException("Family has 0 valid kmers: "+input.family);}
			for(long kmer : seen.toArray()){familyFrequency.increment(kmer);}
		}
		final int workers=Tools.max(1, Tools.min(Shared.threads(), pools.size()));
		final ExecutorService familyExec=Executors.newFixedThreadPool(workers);
		final ArrayList<Future<FamilyResult>> futures=new ArrayList<Future<FamilyResult>>(pools.size());
		for(final FamilyPool pool : pools){
			futures.add(familyExec.submit(new Callable<FamilyResult>(){
				@Override public FamilyResult call(){return selectFamily(pool, familyFrequency);}
			}));
		}
		final ArrayList<FamilyResult> results=new ArrayList<FamilyResult>(pools.size());
		try{
			for(Future<FamilyResult> future : futures){
				try{results.add(future.get());}
				catch(Exception e){throw new RuntimeException("Family covering-set worker failed", e);}
			}
		}finally{familyExec.shutdown();}
		Collections.sort(results, new Comparator<FamilyResult>(){
			@Override public int compare(FamilyResult a, FamilyResult b){return a.family.compareTo(b.family);}
		});
		writeFamilyOutputs(results);
		t.stop(); outstream.println("Processed "+results.size()+" families in one JVM; Time:   "+t);
	}

	private FamilyResult selectFamily(final FamilyPool pool, final LongIntMap familyFrequency){
		final int design=kDesign, use=k;
		final int bits=reducedAlphabet.bits();
		if((long)design*bits>62){throw new RuntimeException("Packed k-mer exceeds 62 bits: k="+design+" bits="+bits+" alphabet="+reducedAlphabet.symbols());}
		final long[][] memberKmers=new long[pool.members.size()][];
		final long[][] memberUniqueKmers=new long[pool.members.size()][];
		final LongIntMap original=new LongIntMap(Math.max(16, pool.members.size()*4));
		int memberIndex=0;
		for(ProteinMember member : pool.members){
			final long[] words=kmers(member.bases, design);
			memberKmers[memberIndex]=words;
			final LongHashSet unique=new LongHashSet(Math.max(16, words.length*2));
			for(long word : words){original.increment(word); unique.add(word);}
			memberUniqueKmers[memberIndex]=unique.toArray();
			memberIndex++;
		}
		if(original.isEmpty()){throw new RuntimeException("Family has 0 valid kmers: "+pool.family);}
		final boolean[] alive=new boolean[pool.members.size()]; Arrays.fill(alive, true);
		final LongHashSet selected=new LongHashSet(Math.max(16, original.size()*2));
		final ArrayList<Selection> selections=new ArrayList<Selection>();
		final LongHashSet rejectedBySpecificity=new LongHashSet(1024);
		int aliveCount=alive.length, rounds=0, currentStep=step;
		while(aliveCount>0 && (maxKmers<=0 || selections.size()<maxKmers)){
			final float coverage=1f-aliveCount/(float)alive.length;
			if(coverage>=minCovFraction){break;}
			final LongIntMap current=new LongIntMap(Math.max(16, original.size()));
			for(int i=0; i<memberKmers.length; i++){
				if(!alive[i]){continue;}
				for(long word : memberKmers[i]){current.increment(word);}
			}
			if(current.isEmpty()){break;}
			final int candidateLimit=Tools.min(current.size(), Tools.min(Integer.MAX_VALUE/2, 2*currentStep));
			final TopKHeap currentHeap=new TopKHeap(candidateLimit);
			final long[] currentKeys=current.keys(); final int[] currentValues=current.values();
			final long currentInvalid=current.invalid();
			for(int cell=0; cell<currentKeys.length; cell++){
				if(currentKeys[cell]!=currentInvalid){currentHeap.add(currentKeys[cell], currentValues[cell], currentKeys[cell]);}
			}
			final long[] candidates=currentHeap.keysDescending();
			final TopKHeap originalHeap=new TopKHeap(Tools.min(currentStep, candidates.length));
			for(long candidate : candidates){originalHeap.add(candidate, original.get(candidate), candidate);}
			final long[] rankedCandidates=originalHeap.keysDescending();
			int added=0;
			for(long word : rankedCandidates){
				if(selected.contains(word)){continue;}
				if(maxFamilies>0 && familyFrequency!=null && familyFrequency.get(word)>maxFamilies){
					rejectedBySpecificity.add(word); continue;
				}
				selected.add(word); selections.add(new Selection(word, original.get(word), rounds+1));
				if(++added>=currentStep || (maxKmers>0 && selections.size()>=maxKmers)){break;}
			}
			if(added==0){break;}
			int evicted=0;
			for(int i=0; i<memberKmers.length; i++){
				if(!alive[i]){continue;}
				int hits=0;
				for(long word : memberUniqueKmers[i]){if(selected.contains(word) && ++hits>=minHits){break;}}
				if(hits>=minHits){alive[i]=false; evicted++;}
			}
			aliveCount-=evicted; rounds++;
				if(evicted==0 && selected.size()>=current.size()){break;}
		}
		final ArrayList<OutputKmer> output=new ArrayList<OutputKmer>();
		final LongHashSet emitted=new LongHashSet(Math.max(16, selections.size()*2));
		int rank=0;
		for(Selection selection : selections){
			final long[] words=useWords(selection.word, design, use);
			for(long word : words){
				if(emitted.add(word)){output.add(new OutputKmer(wordToString(word, use), ++rank, selection.originalCount, selection.round));}
			}
		}
		final int excluded=pool.totalMembers-pool.members.size();
		return new FamilyResult(pool.family, pool.path, pool.inputHash, pool.totalMembers, excluded, output, selections.size(), rounds,
			1f-aliveCount/(float)alive.length, aliveCount, rejectedBySpecificity.size());
	}

	private long[] useWords(final long designWord, final int design, final int use){
		final int bits=reducedAlphabet.bits(), useBits=use*bits;
		final long mask=(1L<<useBits)-1;
		final long[] out=new long[design-use+1];
		for(int i=0; i<out.length; i++){out[i]=(designWord>>(bits*i))&mask;}
		return out;
	}

	private long[] kmers(final byte[] bases, final int length){
		final int bits=reducedAlphabet.bits();
		if((long)length*bits>62){throw new RuntimeException("Packed k-mer exceeds 62 bits: k="+length+" bits="+bits);}
		final long mask=(1L<<(length*bits))-1;
		final LongList out=new LongList(); long word=0; int valid=0;
		for(byte residue : bases){
			final int code=reducedAlphabet.code(residue);
			if(code<0){word=0; valid=0; continue;}
			word=((word<<bits)|code)&mask;
			if(++valid>=length){out.add(word);}
		}
		return out.toArray();
	}

	private String wordToString(long word, final int length){
		final int bits=reducedAlphabet.bits(); final long mask=(1L<<bits)-1;
		final char[] out=new char[length];
		for(int i=length-1; i>=0; i--){out[i]=reducedAlphabet.symbol((int)(word&mask)); word>>>=bits;}
		return new String(out);
	}

	private ArrayList<FamilyInput> loadFamilyInputs(final String path){
		final ArrayList<FamilyInput> out=new ArrayList<FamilyInput>();
		final File f=new File(path);
		if(f.isDirectory()){
			final File[] files=f.listFiles();
			if(files==null){throw new RuntimeException("Cannot read families directory: "+path);}
			Arrays.sort(files, new Comparator<File>(){@Override public int compare(File a, File b){return a.getName().compareTo(b.getName());}});
			for(File child : files){
				if(isProteinFile(child.getName())){
					out.add(new FamilyInput(familyId(child.getName()), child.getPath(), hashFile(child)));
				}
			}
		}else{
			try{
				final BufferedReader br=new BufferedReader(new FileReader(f));
				String line;
				while((line=br.readLine())!=null){
					if(line.length()==0 || line.charAt(0)=='#'){continue;}
					final String[] fields=line.split("\\t", -1);
					if(fields.length!=2){throw new RuntimeException("Malformed family manifest row: "+line);}
					File familyFile=new File(fields[1]);
					if(!familyFile.isAbsolute() && f.getParentFile()!=null){familyFile=new File(f.getParentFile(), fields[1]);}
					out.add(new FamilyInput(fields[0], familyFile.getPath(), hashFile(familyFile)));
				}
				br.close();
			}catch(RuntimeException e){throw e;
			}catch(Exception e){throw new RuntimeException("Could not read family manifest: "+path, e);}
		}
		if(out.isEmpty()){throw new RuntimeException("No family FASTA inputs found: "+path);}
		final HashSet<String> seen=new HashSet<String>();
		for(FamilyInput input : out){
			//A family id is a TSV field and a FASTA-derived name: it may contain '|' (BBTools tid|...| headers, e.g. the mag-qc
			//family rep ids) and any other printable byte; it must not be empty, contain whitespace (tab is the field
			//separator), or start with '#' (a comment line to every consumer). UMP45 2026-09-02: the previous
			//[A-Za-z0-9._-]+ check rejected every real rep id (jobs 25489956-59).
			if(input.family.isEmpty() || input.family.charAt(0)=='#' || input.family.matches(".*\\s.*")){
				throw new RuntimeException("Unsafe family id (empty, leading '#', or whitespace): '"+input.family+"'");
			}
			if(!seen.add(input.family)){throw new RuntimeException("Duplicate family id: "+input.family);}
		}
		return out;
	}

	private static boolean isProteinFile(final String name){
		final String n=name.toLowerCase();
		return n.endsWith(".faa") || n.endsWith(".faa.gz") || n.endsWith(".fa") || n.endsWith(".fa.gz") ||
			n.endsWith(".fasta") || n.endsWith(".fasta.gz");
	}

	private static String familyId(String name){
		final String n=name.replaceFirst("\\.gz$", "");
		return n.replaceFirst("\\.(faa|fa|fasta)$", "");
	}

	private Set<String> loadExclude(final String path){
		final HashSet<String> out=new HashSet<String>();
		if(path==null){return out;}
		try{
			final BufferedReader br=new BufferedReader(new FileReader(path));
			String line;
			while((line=br.readLine())!=null){
				if(line.length()==0 || line.charAt(0)=='#'){continue;}
				final int tab=line.indexOf('\t');
				final String id=(tab<0 ? line : line.substring(0, tab));
				if(id.length()>0 && !out.add(id)){throw new RuntimeException("Duplicate exclude id: "+id);}
			}
			br.close();
		}catch(RuntimeException e){throw e;
		}catch(Exception e){throw new RuntimeException("Could not read exclude file: "+path, e);}
		return out;
	}

	private ProteinLoad loadProteinMembers(final String path, final Set<String> excluded, final String family){
		final File file=new File(path);
		final File source=resolveSource(file);
		final FileFormat ff=FileFormat.testInput(source.getPath(), FileFormat.FASTA, null, true, true);
		final Streamer streamer=StreamerFactory.makeStreamer(ff, 0, true, -1);
		final ArrayList<ProteinMember> out=new ArrayList<ProteinMember>();
		final HashSet<String> ids=new HashSet<String>();
		int total=0;
		streamer.start();
		for(ListNum<Read> ln=streamer.nextList(); ln!=null; ln=streamer.nextList()){
			for(Read read : ln){
				if(read.bases==null || read.bases.length==0){continue;}
				total++;
				final String id=firstToken(read.id);
				if(!ids.add(id)){throw new RuntimeException("Duplicate member id in family "+family+": "+id);}
				Tools.toUpperCase(read.bases);
				if(!excluded.contains(id)){out.add(new ProteinMember(id, read.bases));}
			}
		}
		streamer.close();
		if(streamer.errorState()){throw new RuntimeException("Unreadable family FASTA: "+family+" ["+path+"]");}
		if(total==0){throw new RuntimeException("FASTA contains no non-empty sequences: "+family+" ["+path+"]");}
		return new ProteinLoad(out, total);
	}

	private static String firstToken(final String id){
		int end=id.length();
		for(int i=0; i<id.length(); i++){if(Character.isWhitespace(id.charAt(i))){end=i; break;}}
		return id.substring(0, end);
	}

	private PrintWriter openTextOutput(final String path){
		try{
			final File f=new File(path); final File parent=f.getParentFile();
			if(parent!=null){parent.mkdirs();}
			if(f.exists() && !overwrite){throw new RuntimeException("Output exists and overwrite=f: "+path);}
			return new PrintWriter(new FileOutputStream(f));
		}catch(Exception e){throw new RuntimeException("Could not open output: "+path, e);}
	}

	private void writeFamilyOutputs(final ArrayList<FamilyResult> results){
		final PrintWriter sets=openTextOutput(out);
		sets.println("#schema_version\t1");
		sets.println("#kind\tcovering_set");
		sets.println("#bbtools_commit\t"+bbtoolsCommit);
		sets.println("#alphabet\t"+reducedAlphabet.symbols());
		sets.println("#key\t"+(keySpec==null ? "none" : keySpec));
		sets.println("#k\t"+k+"\tkdesign\t"+kDesign+"\tminhits\t"+minHits+"\tmaxfamilies\t"+maxFamilies);
		sets.println("#columns\tfamily_id\tkmer\tselection_rank\toriginal_count\tround");
		long removed=0;
		for(FamilyResult result : results){
			removed+=result.specificityRemoved;
			for(OutputKmer row : result.output){sets.println(result.family+'\t'+row.text+'\t'+row.rank+'\t'+row.originalCount+'\t'+row.round);}
		}
		sets.println("#cross_family_candidates_removed\t"+removed);
		sets.close();

		final PrintWriter sums=openTextOutput(summary);
		sums.println("#schema_version\t1");
		sums.println("#kind\tcovering_set_summary");
		sums.println("#bbtools_commit\t"+bbtoolsCommit);
		sums.println("#alphabet\t"+reducedAlphabet.symbols());
		sums.println("#key\t"+(keySpec==null ? "none" : keySpec));
		sums.println("#k\t"+k+"\tkdesign\t"+kDesign+"\tminhits\t"+minHits+"\tmaxfamilies\t"+maxFamilies);
		sums.println("#input\tfamily_id\tpath\tsha256");
		for(FamilyResult result : results){sums.println("#input\t"+result.family+'\t'+result.path+'\t'+result.inputHash);}
		sums.println("#cross_family_candidates_removed\t"+removed);
		sums.println("#columns\tfamily_id\tmembers\texcluded\tselected_kmers\trounds\tcoverage\tuncovered");
		for(FamilyResult result : results){
			sums.println(result.family+'\t'+result.members+'\t'+result.excluded+'\t'+result.selectedKmers+'\t'+result.rounds+'\t'+
				String.format(java.util.Locale.ROOT, "%.6f", result.coverage)+'\t'+result.uncovered);
		}
		sums.close();
	}

	private static String hashFile(final File file){
		try{
			final File source=resolveSource(file);
			final MessageDigest digest=MessageDigest.getInstance("SHA-256");
			final FileInputStream in=new FileInputStream(source); final byte[] buffer=new byte[1<<16];
			for(int n=in.read(buffer); n>=0; n=in.read(buffer)){if(n>0){digest.update(buffer, 0, n);}}
			in.close(); final StringBuilder sb=new StringBuilder(64);
			for(byte b : digest.digest()){sb.append(String.format("%02x", b&255));}
			return sb.toString();
		}catch(RuntimeException e){throw e;
		}catch(Exception e){throw new RuntimeException("Could not hash input: "+file, e);}
	}

	/** Resolves a regular file or a symlink chain to its final regular target. */
	private static File resolveSource(final File file){
		try{
			final java.nio.file.Path path=file.toPath();
			if(!java.nio.file.Files.isSymbolicLink(path)){
				if(!java.nio.file.Files.isRegularFile(path)){throw new RuntimeException("Input is not a regular file: "+file);}
				return file;
			}
			final java.nio.file.Path target=path.toRealPath();
			if(!java.nio.file.Files.isRegularFile(target)){
				throw new RuntimeException("Symlink target is not a regular file: "+file+" -> "+target);
			}
			return target.toFile();
		}catch(RuntimeException e){throw e;
		}catch(Exception e){throw new RuntimeException("Unreadable symlink input: "+file, e);}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	private ArrayList<byte[]> loadPool(){
		ArrayList<byte[]> pool=new ArrayList<>();
		loadFasta(in, pool);
		if(extra!=null){
			int base=pool.size();
			ArrayList<byte[]> extraSeqs=new ArrayList<>();
			loadFasta(extra, extraSeqs);
			for(int c=0; c<copies; c++){
				for(byte[] seq : extraSeqs){pool.add(seq);}
			}
			outstream.println("Added "+extraSeqs.size()+" extra sequences x"+copies+
					" copies (pool "+base+" -> "+pool.size()+")");
		}
		return pool;
	}

	/** Loads a FASTA file into memory via Streamer, not ConcurrentReadInputStream.
	 * CRIS spawns exactly one ReadThread for a single unpaired file, so the whole
	 * file's line-splitting/header-parsing/Read-construction runs on ONE core no
	 * matter how many cores are available -- confirmed pinning a single core at
	 * 99% for the entire load phase on a real 2.4Gbp input (Brian, 2026-08-26).
	 * StreamerFactory.makeStreamer(ff, pairnum, ordered, maxReads) with threads
	 * left at its default (-1, "auto") dispatches to the real parallel
	 * FastaStreamer when Shared.threads()>=8: one I/O thread splitting on '>'
	 * boundaries feeds an ordered queue that DEFAULT_THREADS=3 worker threads
	 * drain concurrently, each doing the actual Read construction/validation.
	 * (StreamerFactory.getReads(...) -- the tempting one-liner -- is NOT this:
	 * it hardcodes threads=1, which routes to FastaStreamerST, a single dedicated
	 * worker thread -- same one-core bottleneck as the old CRIS code, just a
	 * newer class name. Call makeStreamer directly to actually get the parallel
	 * path.) Below Shared.threads()<8 this transparently falls back to a
	 * single-threaded streamer, so this is never worse than the old behavior. */
	private void loadFasta(String fname, ArrayList<byte[]> dest){
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTA, null, true, true);
		Streamer st=StreamerFactory.makeStreamer(ff, 0, true, -1);
		st.start();
		for(ListNum<Read> ln=st.nextList(); ln!=null; ln=st.nextList()){
			for(Read r : ln){
				if(r.bases!=null && r.length()>0){
					Tools.toUpperCase(r.bases);
					dest.add(r.bases);
				}
			}
		}
		st.close();
		//Reference-data load: a partial/corrupt read here silently produces a
		//wrong kmer pool with no error visible downstream, so crash loud instead
		//of continuing -- same idiom as StreamerFactory.getReads' own error path.
		if(st.errorState()){
			KillSwitch.kill("Error: a read error (corrupt or truncated input) occurred reading "
				+fname+"; aborting rather than building a covering set from partial/incorrect data.");
		}
	}

	/** Minimum pool size before bothering to fan out across threads -- below this,
	 * per-task submission overhead would exceed the work being parallelized. */
	private static final int MIN_PARALLEL_POOL=64;

	/** Worker pool for the per-round count/evict fan-out, sized by Shared.threads()
	 * and created lazily on first use. REUSED across every round of the greedy
	 * set-cover loop -- spawning fresh Threads per round cost more in thread-startup
	 * overhead than it saved once the algorithm's ~18+ rounds each paid that tax; a
	 * persistent pool amortizes it to (effectively) once per run, matching the
	 * pool()/Future pattern already used for the same reason in
	 * TrnaConsensusBuilder. */
	private ExecutorService pool(){
		if(execPool==null){execPool=Executors.newFixedThreadPool(Tools.max(1, Shared.threads()));}
		return execPool;
	}
	private ExecutorService execPool;

	/** Partitioned kmer counter (Brian's design, 2026-08-26). The previous
	 * multithreaded version gave each thread its own private LongIntMap, then
	 * folded all of them into one map at the end -- correct, but that final
	 * fold is inherently serial and, measured directly, ate roughly a third of
	 * the whole multithreaded runtime (worse the fewer times each kmer repeats,
	 * since repeats are what the fold-by-increment amortizes over). It also caps
	 * total distinct kmers at whatever one LongIntMap can hold.
	 * <p>
	 * This splits the kmer space itself across `numPartitions` independent
	 * LongIntMap shards, each guarded by its own lock. A thread never writes a
	 * shard directly -- it accumulates kmers into its own small per-partition
	 * buffer (a LocalBuffer, NOT shared with other threads) and only takes a
	 * shard's lock to flush a full batch. So lock acquisitions happen once per
	 * `bufferSize` kmers per partition, not once per kmer, keeping contention
	 * low even with many threads hammering the same run. Once every thread's
	 * range is processed and its buffers flushed, the counts are already in
	 * their final home -- there is no merge step at all. Splitting across
	 * shards also means total capacity is N separate LongIntMaps instead of
	 * one, so a corpus with more distinct kmers than a single table comfortably
	 * holds (the motivating case: RefSeq Plastid, 2.4Gbp, ~500M distinct
	 * kmers) is no longer a ceiling. */
	private final class PartitionedCounter {
		final LongIntMap[] shards;
		final Object[] locks;

		PartitionedCounter(long expectedDistinct){
			shards=new LongIntMap[numPartitions];
			locks=new Object[numPartitions];
			final int perShard=(int)Tools.max(256, expectedDistinct/numPartitions+1);
			for(int i=0; i<numPartitions; i++){
				shards[i]=new LongIntMap(perShard);
				locks[i]=new Object();
			}
		}

		private int way(long kmer){return (int)(kmer%numPartitions);}

		/** Per-thread accumulation buffer -- create one per worker via
		 * newLocalBuffer(), never share it across threads, and call flushAll()
		 * exactly once when that worker's assigned range is fully processed. */
		final class LocalBuffer {
			private final long[][] buf=new long[numPartitions][bufferSize];
			private final int[] size=new int[numPartitions];

			void add(long kmer){
				final int w=way(kmer);
				int s=size[w];
				buf[w][s]=kmer;
				s++;
				if(s>=bufferSize){flush(w, s); s=0;}
				size[w]=s;
			}

			private void flush(int w, int n){
				final long[] b=buf[w];
				final LongIntMap shard=shards[w];
				synchronized(locks[w]){
					for(int i=0; i<n; i++){shard.increment(b[i]);}
				}
			}

			void flushAll(){
				for(int w=0; w<numPartitions; w++){
					if(size[w]>0){flush(w, size[w]); size[w]=0;}
				}
			}
		}

		LocalBuffer newLocalBuffer(){return new LocalBuffer();}

		int get(long kmer){return shards[way(kmer)].get(kmer);}

		long size(){
			long sum=0;
			for(LongIntMap m : shards){sum+=m.size();}
			return sum;
		}

		boolean isEmpty(){return size()==0;}
	}

	/** Fills pool[from,to) (or only the alive ones there, if alive!=null) into a
	 * per-thread LocalBuffer. Rolls a forward kmer and (when rcomp is set) a
	 * reverse-complement kmer in lockstep and buffers only the canonical (max
	 * of the two) form -- same idiom as prok/RiboMaker.loadFilter -- so a motif
	 * and its reverse complement are counted as the SAME kmer instead of two
	 * independent ones. `buf` must be this thread's own LocalBuffer: safe to
	 * call from multiple threads concurrently as long as each call gets a
	 * distinct buffer and a disjoint [from,to) range of the same pool/alive. */
	private void fillRange(ArrayList<byte[]> pool, boolean[] alive, int from, int to,
			int kLen, long mask, PartitionedCounter.LocalBuffer buf){
		final int shift2=2*kLen-2;
		for(int i=from; i<to; i++){
			if(alive!=null && !alive[i]){continue;}
			byte[] bases=pool.get(i);
			long kmer=0, rkmer=0;
			int len=0;
			for(byte b : bases){
				int num=AminoAcid.baseToNumber[b];
				if(num>=0){
					kmer=((kmer<<2)|num)&mask;
					if(rcomp){
						long comp=AminoAcid.baseToComplementNumber[b];
						rkmer=((rkmer>>>2)|(comp<<shift2))&mask;
					}
					len++;
					if(len>=kLen){buf.add(rcomp ? Tools.max(kmer, rkmer) : kmer);}
				}else{len=0; kmer=0; rkmer=0;}
			}
		}
		buf.flushAll();
	}

	/** Multithreaded partitioned kmer counting -- see PartitionedCounter's own
	 * doc for why this replaced the previous private-map-per-thread-then-merge
	 * design. alive==null counts the whole pool (the one-time full-corpus
	 * count); a non-null alive[] counts only the currently-uncovered sequences
	 * (the per-round recount). expectedDistinct is only a SIZING HINT for the
	 * shards (never a correctness constraint) -- pass the best upper-bound
	 * estimate available; under- or over-estimating just costs a few extra
	 * resizes or a little unused capacity. Falls back to a single-threaded fill
	 * (still partitioned, so the return type and downstream code are uniform)
	 * when the pool is too small to be worth splitting across threads. */
	private PartitionedCounter countKmersMT(final ArrayList<byte[]> pool, final boolean[] alive,
			final int kLen, final long mask, long expectedDistinct){
		final PartitionedCounter counter=new PartitionedCounter(expectedDistinct);
		final int n=pool.size();
		final int threads=Tools.max(1, Shared.threads());
		if(threads<=1 || n<MIN_PARALLEL_POOL){
			fillRange(pool, alive, 0, n, kLen, mask, counter.newLocalBuffer());
			return counter;
		}
		final int chunk=(n+threads-1)/threads;
		final ArrayList<Future<?>> futures=new ArrayList<>(threads);
		for(int t=0; t<threads; t++){
			final int from=t*chunk, to=Tools.min(n, from+chunk);
			if(from>=to){continue;}
			final PartitionedCounter.LocalBuffer buf=counter.newLocalBuffer();
			futures.add(pool().submit(new Runnable(){
				@Override
				public void run(){fillRange(pool, alive, from, to, kLen, mask, buf);}
			}));
		}
		for(Future<?> f : futures){
			try{f.get();}catch(Exception e){throw new RuntimeException(e);}
		}
		return counter;
	}

	/** Returns the top N kmers by count without materializing or sorting the
	 * complete kmer space. Each shard is scanned independently in parallel into
	 * a bounded primitive heap, then the small shard results are merged. Numeric
	 * kmer order breaks count ties reproducibly, independent of shard count. */
	private long[] topNByCount(final PartitionedCounter counter, final int n){
		if(n<=0 || counter.isEmpty()){return new long[0];}
		final int limit=(int)Tools.min(n, counter.size());
		final int ways=counter.shards.length;
		final int threads=Tools.max(1, Shared.threads());
		final TopKHeap merged=new TopKHeap(limit);
		if(threads<=1 || ways<=1){
			for(int way=0; way<ways; way++){merged.add(topKFromShard(counter.shards[way], limit));}
		}else{
			final ArrayList<Future<TopKHeap>> futures=new ArrayList<>(ways);
			for(int way=0; way<ways; way++){
				final int w=way;
				futures.add(pool().submit(new Callable<TopKHeap>(){
					@Override
					public TopKHeap call(){return topKFromShard(counter.shards[w], limit);}
				}));
			}
			for(Future<TopKHeap> f : futures){
				try{merged.add(f.get());}catch(Exception e){throw new RuntimeException(e);}
			}
		}
		return merged.keysDescending();
	}

	/** Scans one immutable shard into a bounded top-K heap. */
	private static TopKHeap topKFromShard(LongIntMap map, int limit){
		final TopKHeap heap=new TopKHeap(Tools.min(limit, map.size()));
		final long[] keys=map.keys();
		final int[] counts=map.values();
		final long invalid=map.invalid();
		for(int cell=0; cell<keys.length; cell++){
			final long key=keys[cell];
			if(key!=invalid){heap.add(key, counts[cell], key);}
		}
		return heap;
	}

	/** From the candidate list, re-rank by original prevalence, keep top step. */
	private long[] rankByOriginal(long[] candidates, PartitionedCounter originalCounts, int keep){
		final TopKHeap heap=new TopKHeap(Tools.min(keep, candidates.length));
		for(int i=0; i<candidates.length; i++){
			heap.add(candidates[i], originalCounts.get(candidates[i]), candidates[i]);
		}
		return heap.keysDescending();
	}

	/** Projects the next batch size from the previous round's observed
	 * sequences-per-kmer recovery. Bounds and the optional round-2 boost are
	 * tunable because conservation changes the useful step scale. */
	private int adaptiveStep(int previousKmers, int previousRecovered, int remaining, int nextRound){
		double projected=previousKmers*(stepFraction*remaining)/Tools.max(1, previousRecovered);
		if(nextRound==2){projected*=step2Boost;}
		final long raw=projected>=Integer.MAX_VALUE ? Integer.MAX_VALUE : (long)Math.ceil(projected);
		final long min=scaledStep(minStepMult);
		final long max=Tools.max(min, scaledStep(maxStepMult));
		return (int)Tools.mid(min, max, raw);
	}

	private long scaledStep(double mult){
		final double x=Math.ceil(step*mult);
		return Tools.max(1L, x>=Integer.MAX_VALUE ? Integer.MAX_VALUE : (long)x);
	}

	/** Marks sequences as not-alive if they contain any selected (canonical) kmer,
	 * over pool[from,to). Only reads `selected` (a fixed snapshot for the whole
	 * round) and writes exclusively to its own index range of `alive`, so this is
	 * safe to run from multiple threads on disjoint ranges concurrently.
	 * Returns the number of newly evicted sequences in this range. */
	private int evictRange(ArrayList<byte[]> pool, boolean[] alive,
			LongHashSet selected, int from, int to, int kLen, long mask){
		int evicted=0;
		final int shift2=2*kLen-2;
		for(int i=from; i<to; i++){
			if(!alive[i]){continue;}
			byte[] bases=pool.get(i);
			long kmer=0, rkmer=0;
			int len=0;
			for(byte b : bases){
				int num=AminoAcid.baseToNumber[b];
				if(num>=0){
					kmer=((kmer<<2)|num)&mask;
					if(rcomp){
						long comp=AminoAcid.baseToComplementNumber[b];
						rkmer=((rkmer>>>2)|(comp<<shift2))&mask;
					}
					len++;
					if(len>=kLen && selected.contains(rcomp ? Tools.max(kmer, rkmer) : kmer)){
						alive[i]=false;
						evicted++;
						break;
					}
				}else{len=0; kmer=0; rkmer=0;}
			}
		}
		return evicted;
	}

	/** Multithreaded wrapper for evictRange -- same chunking/fallback/reused-pool
	 * strategy as countKmersMT. Threads write disjoint indices of `alive`, so no
	 * merge step is needed beyond summing the per-thread eviction counts. */
	private int evictCoveredMT(final ArrayList<byte[]> pool, final boolean[] alive,
			final LongHashSet selected, final int kLen, final long mask){
		final int n=pool.size();
		final int threads=Tools.max(1, Shared.threads());
		if(threads<=1 || n<MIN_PARALLEL_POOL){
			return evictRange(pool, alive, selected, 0, n, kLen, mask);
		}
		final int chunk=(n+threads-1)/threads;
		final ArrayList<Future<Integer>> futures=new ArrayList<>(threads);
		for(int t=0; t<threads; t++){
			final int from=t*chunk, to=Tools.min(n, from+chunk);
			if(from>=to){continue;}
			futures.add(pool().submit(new Callable<Integer>(){
				@Override
				public Integer call(){return evictRange(pool, alive, selected, from, to, kLen, mask);}
			}));
		}
		int total=0;
		for(Future<Integer> f : futures){
			try{total+=f.get();}catch(Exception e){throw new RuntimeException(e);}
		}
		return total;
	}

	/** Converts selected kmer longs into single-kmer FASTA reads. */
	private ArrayList<Read> kmersToReads(LongList kmers, int kLen){
		ArrayList<Read> list=new ArrayList<>(kmers.size);
		for(int i=0; i<kmers.size; i++){
			byte[] bases=longToBytes(kmers.get(i), kLen);
			list.add(new Read(bases, null, "kmer_"+i, i));
		}
		return list;
	}

	/** For k+1 design: each selected (k+1)-mer yields two overlapping k-mers.
	 * Deduplicate and output only the unique k-mers, preserving selection order. */
	private ArrayList<Read> designToUseKmers(LongList designKmers, int kD, int kU){
		final long useMask=~((-1L)<<(2*kU));
		LongHashSet seen=new LongHashSet(designKmers.size*2);
		LongList useKmers=new LongList();
		for(int i=0; i<designKmers.size; i++){
			long dkmer=designKmers.get(i);
			// A (kD)-mer contains (kD-kU+1) overlapping (kU)-mers
			for(int offset=0; offset<=kD-kU; offset++){
				long ukmer=(dkmer>>(2*offset))&useMask;
				if(!seen.contains(ukmer)){
					seen.add(ukmer);
					useKmers.add(ukmer);
				}
			}
		}
		outstream.println("Design "+kD+"-mers: "+designKmers.size+" -> use "+kU+"-mers: "+useKmers.size);
		ArrayList<Read> list=new ArrayList<>(useKmers.size);
		for(int i=0; i<useKmers.size; i++){
			byte[] bases=longToBytes(useKmers.get(i), kU);
			list.add(new Read(bases, null, "kmer_"+i, i));
		}
		return list;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Utility Methods       ----------------*/
	/*--------------------------------------------------------------*/

	private static byte[] longToBytes(long kmer, int kLen){
		byte[] bases=new byte[kLen];
		for(int i=kLen-1; i>=0; i--){
			bases[i]=AminoAcid.numberToBase[(int)(kmer&3)];
			kmer>>>=2;
		}
		return bases;
	}

	/** Primitive bounded min-heap. The root is the worst retained item: lower
	 * count, then later stable order. This avoids boxed entries in the hot rank. */
	private static final class TopKHeap {
		TopKHeap(int capacity_){
			capacity=capacity_;
			keys=new long[capacity];
			counts=new int[capacity];
			orders=new long[capacity];
		}

		void add(TopKHeap b){
			for(int i=0; i<b.size; i++){add(b.keys[i], b.counts[i], b.orders[i]);}
		}

		boolean add(long key, int count, long order){
			if(capacity<=0){return false;}
			if(size<capacity){
				keys[size]=key;
				counts[size]=count;
				orders[size]=order;
				siftUp(size);
				size++;
				return true;
			}
			if(!better(count, order, counts[0], orders[0])){return false;}
			keys[0]=key;
			counts[0]=count;
			orders[0]=order;
			siftDown(0);
			return true;
		}

		long[] keysDescending(){
			final long[] out=new long[size];
			for(int i=size-1; i>=0; i--){
				out[i]=keys[0];
				size--;
				if(size>0){
					keys[0]=keys[size];
					counts[0]=counts[size];
					orders[0]=orders[size];
					siftDown(0);
				}
			}
			return out;
		}

		private void siftUp(int child){
			final long key=keys[child], order=orders[child];
			final int count=counts[child];
			while(child>0){
				final int parent=(child-1)>>>1;
				if(!worse(count, order, counts[parent], orders[parent])){break;}
				keys[child]=keys[parent];
				counts[child]=counts[parent];
				orders[child]=orders[parent];
				child=parent;
			}
			keys[child]=key;
			counts[child]=count;
			orders[child]=order;
		}

		private void siftDown(int parent){
			final long key=keys[parent], order=orders[parent];
			final int count=counts[parent];
			for(int child=parent*2+1; child<size; child=parent*2+1){
				if(child+1<size && worse(counts[child+1], orders[child+1],
						counts[child], orders[child])){child++;}
				if(!worse(counts[child], orders[child], count, order)){break;}
				keys[parent]=keys[child];
				counts[parent]=counts[child];
				orders[parent]=orders[child];
				parent=child;
			}
			keys[parent]=key;
			counts[parent]=count;
			orders[parent]=order;
		}

		private static boolean better(int ac, long ao, int bc, long bo){
			return ac>bc || (ac==bc && ao<bo);
		}

		private static boolean worse(int ac, long ao, int bc, long bo){
			return ac<bc || (ac==bc && ao>bo);
		}

		private final int capacity;
		private final long[] keys;
		private final int[] counts;
		private final long[] orders;
		private int size=0;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String in;
	private String out;
	private String extra;
	private String families;
	private String exclude;
	private String summary;
	private String alphabetSpec;
	private String keySpec;
	private String bbtoolsCommit="unknown";
	private ReducedAlphabet reducedAlphabet;
	private boolean proteinMode=false;
	private final FileFormat ffin;
	private final FileFormat ffout;
	private boolean overwrite=true;
	private boolean errorState=false;

	private int k=17;
	private int kDesign=-1;
	private int step=500;
	/** Fraction of remaining sequences targeted per adaptive round; 0 keeps the
	 * historical fixed-step behavior. */
	private double stepFraction=0;
	private double minStepMult=0.1;
	private double maxStepMult=2;
	private double step2Boost=1;
	private int maxKmers=0;
	private int minHits=1;
	private int maxFamilies=0;
	private float minCovFraction=0.999f;
	private int copies=10;
	private boolean rcomp=false;
	/** Number of PartitionedCounter shards. <=0 (default) resolves in process()
	 * to the nearest odd value >= max(15, 2*Shared.threads()) -- Brian's
	 * illustrative floor of 15, doubled against thread count so flush-time lock
	 * collisions between threads stay rare, and forced odd so kmer%numPartitions
	 * doesn't alias against 2-bit-packed kmers' structured low bits. */
	private int numPartitions=-1;
	/** Per-thread, per-partition batch size before a flush takes that
	 * partition's lock (Brian's illustrative 200). */
	private int bufferSize=200;

	private PrintStream outstream=System.err;
	private static boolean verbose=false;

	private static final class FamilyInput {
		FamilyInput(String family_, String path_, String sha256_){family=family_; path=path_; sha256=sha256_;}
		final String family, path, sha256;
	}
	private static final class ProteinMember {
		ProteinMember(String id_, byte[] bases_){id=id_; bases=bases_;}
		final String id; final byte[] bases;
	}
	private static final class ProteinLoad {
		ProteinLoad(ArrayList<ProteinMember> members_, int total_){members=members_; totalMembers=total_;}
		final ArrayList<ProteinMember> members; final int totalMembers;
	}
	private static final class FamilyPool {
		FamilyPool(String family_, ArrayList<ProteinMember> members_, String path_, String hash_, int total_){
			family=family_; members=members_; path=path_; inputHash=hash_; totalMembers=total_;
		}
		final String family, path, inputHash; final ArrayList<ProteinMember> members; final int totalMembers;
	}
	private static final class Selection {
		Selection(long word_, int originalCount_, int round_){word=word_; originalCount=originalCount_; round=round_;}
		final long word; final int originalCount, round;
	}
	private static final class OutputKmer {
		OutputKmer(String text_, int rank_, int originalCount_, int round_){text=text_; rank=rank_; originalCount=originalCount_; round=round_;}
		final String text; final int rank, originalCount, round;
	}
	private static final class FamilyResult {
		FamilyResult(String family_, String path_, String hash_, int members_, int excluded_, ArrayList<OutputKmer> output_,
				int selectedKmers_, int rounds_, float coverage_, int uncovered_, int specificityRemoved_){
			family=family_; path=path_; inputHash=hash_; members=members_; excluded=excluded_; output=output_;
			selectedKmers=selectedKmers_; rounds=rounds_; coverage=coverage_; uncovered=uncovered_; specificityRemoved=specificityRemoved_;
		}
		final String family, path, inputHash; final int members, excluded, selectedKmers, rounds, uncovered, specificityRemoved;
		final float coverage; final ArrayList<OutputKmer> output;
	}
}
