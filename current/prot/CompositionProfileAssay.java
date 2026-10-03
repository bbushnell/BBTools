package prot;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import structures.ByteBuilder;

/**
 * Composition-profile shortlist assay (Brian, 2026-09-03): builds one amino-acid k-mer
 * composition profile per protein family (k=1: 20 or 14 dims; k=2: 400 or 196 dims,
 * over the aa20 or c14 alphabet), then scores every sampled query against every family
 * by cosine similarity to see whether composition alone is a usable shortlist signal
 * (cheap relative to seeding+alignment).
 *
 * <p>This is an evidence run, not production code (spec: "answer first, engineer
 * later"). A family's profile is the SUM of its members' residue (k=1) or adjacent-pair
 * (k=2) counts; the held-out profile set subtracts each sampled query's own counts from
 * its family's total, so scoring never sees the query inside its own profile. The
 * background is the corpus-wide k-mer frequency (computed once, shared by both profile
 * sets and every variant); {@code background=members} weights it by member count
 * (a few huge families can dominate it), {@code background=families} weights every
 * family equally (mean of the 4,432 per-family frequency vectors).</p>
 *
 * <p>Smoothing is variant-specific (UMP45, 2026-09-03): cosine similarity on raw counts
 * needs no smoothing, and a pseudocount comparable to or larger than a short query's
 * total count (e.g. 1.0 with 400 dims against a ~100-residue k=2 query) swamps the real
 * signal with flat pseudocount mass. {@code rawpseudo=} (default 0) feeds the {@code raw}
 * and {@code centered} variants; {@code logpseudo=} (default 1e-3, must be positive) feeds
 * only {@code logodds}, where a pseudocount is required to keep {@code log(freq/bgFreq)}
 * finite.</p>
 *
 * <p>Usage: {@code java -ea prot.CompositionProfileAssay familiesdir=<dir>
 * ranktorepid=<rank_to_repid.tsv> queries=<queries.faa> truth=<truth.tsv> out=<summary.tsv>
 * [perquery=<perquery.tsv>] [strata=<strata.tsv>] [alphabets=aa20,c14] [ks=1,2]
 * [profilesets=heldout,insample] [variants=raw,centered,logodds] [ns=1,5,25,50,100,1000]
 * [rawpseudo=0] [logpseudo=1e-3] [background=members|families] [maxfamilies=N]}</p>
 *
 * @author Ady
 */
public final class CompositionProfileAssay {

	public static void main(final String[] args){new CompositionProfileAssay(args).process();}

	CompositionProfileAssay(final String[] args){
		for(String arg : args){
			final int eq=arg.indexOf('=');
			final String a=(eq<0 ? arg : arg.substring(0, eq)).toLowerCase();
			final String b=(eq<0 ? null : arg.substring(eq+1));
			if(a.equals("familiesdir") || a.equals("families")){familiesDir=b;}
			else if(a.equals("ranktorepid") || a.equals("rank2repid")){rankToRepIdFile=b;}
			else if(a.equals("queries") || a.equals("query")){queryFile=b;}
			else if(a.equals("truth")){truthFile=b;}
			else if(a.equals("out")){outFile=b;}
			else if(a.equals("perquery")){perQueryFile=b;}
			else if(a.equals("strata")){strataFile=b;}
			else if(a.equals("alphabets")){alphabetNames=split(b);}
			else if(a.equals("ks")){ks=parseInts(b);}
			else if(a.equals("profilesets")){profileSetNames=split(b);}
			else if(a.equals("variants")){variantNames=split(b);}
			else if(a.equals("ns") || a.equals("topn")){topNs=parseInts(b);}
			else if(a.equals("rawpseudo")){rawPseudocount=Double.parseDouble(b);}
			else if(a.equals("logpseudo")){logPseudocount=Double.parseDouble(b);}
			else if(a.equals("background")){backgroundMode=b.toLowerCase();}
			else if(a.equals("maxfamilies")){maxFamilies=Integer.parseInt(b);}
			else{throw new RuntimeException("Unknown argument: "+arg);}
		}
		if(familiesDir==null || rankToRepIdFile==null || queryFile==null || truthFile==null || outFile==null){
			throw new RuntimeException("familiesdir=, ranktorepid=, queries=, truth=, and out= are required.");
		}
		Arrays.sort(topNs);
		for(int n : topNs){if(n<1){throw new RuntimeException("top-N must be positive: "+n);}}
		for(int k : ks){if(k!=1 && k!=2){throw new RuntimeException("Only k=1 and k=2 are implemented: "+k);}}
		if(rawPseudocount<0){throw new RuntimeException("rawpseudo must be >=0: "+rawPseudocount);}
		if(logPseudocount<=0){throw new RuntimeException("logpseudo must be positive (used to avoid log(0)): "+logPseudocount);}
		if(!backgroundMode.equals("members") && !backgroundMode.equals("families")){
			throw new RuntimeException("background must be 'members' or 'families': "+backgroundMode);
		}
	}

	private void process(){
		final String[] repIds=loadRankToRepId(rankToRepIdFile, maxFamilies);
		final int nFamilies=repIds.length;
		final long[] familyHash=new long[nFamilies];
		final HashMap<String, Integer> familyIndex=new HashMap<String, Integer>(nFamilies*2);
		for(int i=0; i<nFamilies; i++){familyIndex.put(repIds[i], i); familyHash[i]=ReducedAlphabetSeedAssay.stableHash(repIds[i]);}

		final HashMap<String, String> truth=ReducedAlphabetSeedAssay.loadMap(truthFile, "query truth");
		final ArrayList<ProteinSequence> allQueries=ReducedAlphabetSeedAssay.loadFasta(queryFile);
		//Restrict to queries whose truth family survived maxfamilies truncation (precheck mode only; full runs keep all).
		final ArrayList<ProteinSequence> queries=new ArrayList<ProteinSequence>(allQueries.size());
		final int[] queryFamily0=new int[allQueries.size()];
		int nq=0;
		final java.util.HashSet<String> seenIds=new java.util.HashSet<String>();
		for(ProteinSequence q : allQueries){
			if(!seenIds.add(q.id)){throw new RuntimeException("Duplicate query id: "+q.id);}
			final String family=truth.get(q.id);
			if(family==null){throw new RuntimeException("Missing truth for query: "+q.id);}
			final Integer idx=familyIndex.get(family);
			if(idx==null){
				if(maxFamilies>0){continue;}//precheck subset: family was truncated away
				throw new RuntimeException("Truth family has no rank entry: query="+q.id+", family="+family);
			}
			queries.add(q); queryFamily0[nq]=idx; nq++;
		}
		final int[] queryFamily=Arrays.copyOf(queryFamily0, nq);
		System.err.println("Loaded "+nFamilies+" families, "+queries.size()+" queries (of "+allQueries.size()+" sampled).");

		final AlphaK[] configs=buildConfigs(alphabetNames, ks);

		final FileFormat ffOut=FileFormat.testOutput(outFile, FileFormat.TXT, null, true, true, false, false);
		final ByteStreamWriter bsw=new ByteStreamWriter(ffOut);
		bsw.start();
		final ByteStreamWriter pqw=(perQueryFile==null ? null :
			new ByteStreamWriter(FileFormat.testOutput(perQueryFile, FileFormat.TXT, null, true, true, false, false)));
		final ByteStreamWriter stw=(strataFile==null ? null :
			new ByteStreamWriter(FileFormat.testOutput(strataFile, FileFormat.TXT, null, true, true, false, false)));
		if(pqw!=null){pqw.start();}
		if(stw!=null){stw.start();}
		try{
			final ByteBuilder bb=new ByteBuilder(4096);
			appendSummaryHeader(bb); bsw.print(bb); bb.clear();
			if(pqw!=null){appendPerQueryHeader(bb); pqw.print(bb); bb.clear();}
			if(stw!=null){appendStrataHeader(bb); stw.print(bb); bb.clear();}

			for(AlphaK config : configs){
				System.err.println("Scanning family corpus for "+config.label+" (dims="+config.dims+")...");
				final long scanStart=System.nanoTime();
				final double[][] famCounts=new double[nFamilies][config.dims];
				final int[] famMemberCount=new int[nFamilies];
				scanFamilyCorpus(familiesDir, nFamilies, config, famCounts, famMemberCount);
				final long scanNanos=System.nanoTime()-scanStart;
				System.err.println("Scanned in "+(scanNanos/1e9)+"s.");

				final double[][] qCounts=new double[queries.size()][config.dims];
				for(int i=0; i<queries.size(); i++){addCounts(queries.get(i).enc, config.alphabet, config.k, qCounts[i]);}

				final double[][] queryContribByFamily=new double[nFamilies][config.dims];
				final int[] queryCountByFamily=new int[nFamilies];
				for(int i=0; i<queries.size(); i++){
					final double[] dst=queryContribByFamily[queryFamily[i]];
					final double[] src=qCounts[i];
					for(int d=0; d<config.dims; d++){dst[d]+=src[d];}
					queryCountByFamily[queryFamily[i]]++;
				}
				//Held-out rows report the ACTUAL held-out member count (UMP45, 2026-09-03), not the in-sample
				//count -- famMemberCount alone would silently overstate family_size for held-out profileset rows.
				final int[] heldOutMemberCount=new int[nFamilies];
				for(int r=0; r<nFamilies; r++){heldOutMemberCount[r]=famMemberCount[r]-queryCountByFamily[r];}

				//Background: computed once from the full (in-sample) corpus, shared by every profileset/variant
				//(spec: "compute once over all members"). 'members' weights by member count -- for a family that
				//is a large share of the corpus this makes fam/bg~1 for that family's OWN logodds, i.e. near-zero
				//signal (UMP45's diagnosis of the large-family logodds collapse). 'families' weights every family
				//equally (mean of the 4,432 per-family frequency vectors), avoiding that self-referential bias.
				final double[] bgFreqRaw=new double[config.dims], bgFreqLog=new double[config.dims];
				if(backgroundMode.equals("members")){
					final double[] bgCounts=new double[config.dims];
					for(int r=0; r<nFamilies; r++){
						for(int d=0; d<config.dims; d++){bgCounts[d]+=famCounts[r][d];}
					}
					System.arraycopy(toFreq(bgCounts, config.dims, rawPseudocount), 0, bgFreqRaw, 0, config.dims);
					System.arraycopy(toFreq(bgCounts, config.dims, logPseudocount), 0, bgFreqLog, 0, config.dims);
				}else{
					for(int r=0; r<nFamilies; r++){
						final double[] fr=toFreq(famCounts[r], config.dims, rawPseudocount);
						final double[] fl=toFreq(famCounts[r], config.dims, logPseudocount);
						for(int d=0; d<config.dims; d++){bgFreqRaw[d]+=fr[d]; bgFreqLog[d]+=fl[d];}
					}
					for(int d=0; d<config.dims; d++){bgFreqRaw[d]/=nFamilies; bgFreqLog[d]/=nFamilies;}
				}

				for(String profileSetName : profileSetNames){
					final boolean heldOut=profileSetName.equals("heldout");
					if(!heldOut && !profileSetName.equals("insample")){
						throw new RuntimeException("Unknown profileset: "+profileSetName);
					}
					final double[][] sourceCounts=new double[nFamilies][config.dims];
					for(int r=0; r<nFamilies; r++){
						for(int d=0; d<config.dims; d++){
							sourceCounts[r][d]=famCounts[r][d]-(heldOut ? queryContribByFamily[r][d] : 0);
							//Runtime tripwire (UMP45, 2026-09-03): a query's own composition can never exceed its
							//family's total unless the held-out subtraction assumption breaks (queries.faa content
							//diverging from the family .faa copy it was sampled from, or a double-counted query).
							//A negative bin means that assumption is violated -- crash loud rather than score on
							//silently corrupted held-out profiles.
							assert(sourceCounts[r][d]>=-1e-6) : "Held-out count went negative for family rank "+r+
								", dim "+d+", config "+config.label+": famCount="+famCounts[r][d]+" queryContrib="+
								queryContribByFamily[r][d]+" (query composition exceeds its family's total -- the "+
								"id-based held-out removal assumes queries.faa content matches the family .faa copy).";
						}
					}
					final double[][] famFreqRaw=new double[nFamilies][], famFreqLog=new double[nFamilies][];
					for(int r=0; r<nFamilies; r++){
						famFreqRaw[r]=toFreq(sourceCounts[r], config.dims, rawPseudocount);
						famFreqLog[r]=toFreq(sourceCounts[r], config.dims, logPseudocount);
					}
					final double[][] qFreqRaw=new double[queries.size()][], qFreqLog=new double[queries.size()][];
					for(int i=0; i<queries.size(); i++){
						qFreqRaw[i]=toFreq(qCounts[i], config.dims, rawPseudocount);
						qFreqLog[i]=toFreq(qCounts[i], config.dims, logPseudocount);
					}
					final int[] familySize=(heldOut ? heldOutMemberCount : famMemberCount);

					for(String variantName : variantNames){
						final long scoreStart=System.nanoTime();
						final boolean isLog=variantName.equals("logodds");
						final double[][] srcFamFreq=(isLog ? famFreqLog : famFreqRaw);
						final double[][] srcQFreq=(isLog ? qFreqLog : qFreqRaw);
						final double[] srcBg=(isLog ? bgFreqLog : bgFreqRaw);
						final double[][] famVec=new double[nFamilies][config.dims];
						for(int r=0; r<nFamilies; r++){applyVariant(srcFamFreq[r], srcBg, variantName, famVec[r]); normalize(famVec[r]);}
						final double[][] qVec=new double[queries.size()][config.dims];
						for(int i=0; i<queries.size(); i++){applyVariant(srcQFreq[i], srcBg, variantName, qVec[i]); normalize(qVec[i]);}

						final Result result=score(qVec, famVec, queryFamily, familyHash, repIds, config, profileSetName,
							variantName, topNs, queries, familySize, pqw);
						result.append(bb, topNs); bsw.print(bb); bb.clear();
						if(stw!=null){
							appendStrata(stw, bb, result, config, profileSetName, variantName, topNs);
						}
						System.err.println(config.label+" "+profileSetName+" "+variantName+": scored in "+
							((System.nanoTime()-scoreStart)/1e9)+"s, recall@"+topNs[0]+"="+
							(result.recall[0]/(double)result.queries));
					}
				}
			}
		}finally{
			bsw.poisonAndWait();
			if(pqw!=null){pqw.poisonAndWait();}
			if(stw!=null){stw.poisonAndWait();}
		}
	}

	/** One label per (alphabet, k) combination requested; drives dims and the family-corpus scan.
	 *  Package-private (not private) so {@link HybridShortlistAssay} can reuse the corpus scan
	 *  without re-implementing it (2026-09-03). */
	static final class AlphaK {
		AlphaK(final String alphabetName_, final ReducedAlphabet alphabet_, final int k_){
			alphabetName=alphabetName_; alphabet=alphabet_; k=k_;
			dims=(int)Math.round(Math.pow(alphabet.classes(), k));
			label=alphabetName+"_k"+k;
		}
		final String alphabetName, label;
		final ReducedAlphabet alphabet;
		final int k, dims;
	}

	private static AlphaK[] buildConfigs(final String[] alphabetNames, final int[] ks){
		final ArrayList<AlphaK> list=new ArrayList<AlphaK>();
		for(String name : alphabetNames){
			final ReducedAlphabet alphabet=ReducedAlphabet.named(name);
			for(int k : ks){list.add(new AlphaK(name, alphabet, k));}
		}
		return list.toArray(new AlphaK[0]);
	}

	/** Accumulates k=1 residue counts or k=2 adjacent-pair counts (no windowing across sequence boundaries). */
	static void addCounts(final byte[] enc, final ReducedAlphabet alphabet, final int k, final double[] counts){
		final int classes=alphabet.classes();
		if(k==1){
			for(final byte b : enc){
				final int code=alphabet.codeEncoded(b);
				if(code>=0){counts[code]++;}
			}
		}else if(k==2){
			int prev=-1;
			for(final byte b : enc){
				final int code=alphabet.codeEncoded(b);
				if(prev>=0 && code>=0){counts[prev*classes+code]++;}
				prev=code;
			}
		}else{
			throw new RuntimeException("Unsupported k: "+k);
		}
	}

	/** freq[d] = (counts[d]+pseudocount) / (total+dims*pseudocount); always strictly positive. */
	static double[] toFreq(final double[] counts, final int dims, final double pseudocount){
		double total=0;
		for(double c : counts){total+=c;}
		final double denom=total+dims*pseudocount;
		final double[] freq=new double[dims];
		for(int d=0; d<dims; d++){freq[d]=(counts[d]+pseudocount)/denom;}
		return freq;
	}

	/** Background frequency vector for one pseudocount (2026-09-03, pure addition for {@link HybridShortlistAssay}
	 *  reuse -- {@code process()}'s inline background block is untouched to avoid any risk to the gated v3
	 *  pipeline). {@code background="members"} weights by member count (sum of all family counts, then one
	 *  toFreq); {@code "families"} weights every family equally (mean of each family's own frequency vector) --
	 *  see the class javadoc for why "families" avoids a large family dominating its own background. */
	static double[] computeBackgroundFreq(final double[][] famCounts, final int nFamilies, final int dims,
			final double pseudocount, final String backgroundMode){
		final double[] bgFreq=new double[dims];
		if(backgroundMode.equals("members")){
			final double[] bgCounts=new double[dims];
			for(int r=0; r<nFamilies; r++){
				for(int d=0; d<dims; d++){bgCounts[d]+=famCounts[r][d];}
			}
			System.arraycopy(toFreq(bgCounts, dims, pseudocount), 0, bgFreq, 0, dims);
		}else if(backgroundMode.equals("families")){
			for(int r=0; r<nFamilies; r++){
				final double[] fr=toFreq(famCounts[r], dims, pseudocount);
				for(int d=0; d<dims; d++){bgFreq[d]+=fr[d];}
			}
			for(int d=0; d<dims; d++){bgFreq[d]/=nFamilies;}
		}else{
			throw new RuntimeException("background must be 'members' or 'families': "+backgroundMode);
		}
		return bgFreq;
	}

	/** Writes the variant-transformed vector into {@code dst} (raw / centered / logodds vs the shared background). */
	static void applyVariant(final double[] freq, final double[] bgFreq, final String variant, final double[] dst){
		if(variant.equals("raw")){
			System.arraycopy(freq, 0, dst, 0, freq.length);
		}else if(variant.equals("centered")){
			for(int d=0; d<freq.length; d++){dst[d]=freq[d]-bgFreq[d];}
		}else if(variant.equals("logodds")){
			for(int d=0; d<freq.length; d++){dst[d]=Math.log(freq[d]/bgFreq[d]);}
		}else{
			throw new RuntimeException("Unknown variant: "+variant);
		}
	}

	/** L2-normalizes in place; a (pathological) zero-norm vector is left as all-zero, scoring cosine 0 with everything. */
	static void normalize(final double[] v){
		double ss=0;
		for(double x : v){ss+=x*x;}
		final double norm=Math.sqrt(ss);
		if(norm>1e-30){for(int d=0; d<v.length; d++){v[d]/=norm;}}
	}

	/** Package-private + static (2026-09-03): takes every input as a parameter already, so it never needed
	 *  instance state; widened so {@link HybridShortlistAssay} can reuse it instead of re-scanning the corpus. */
	static void scanFamilyCorpus(final String familiesDir, final int nFamilies, final AlphaK config,
			final double[][] famCounts, final int[] famMemberCount){
		for(int rank=0; rank<nFamilies; rank++){
			final String file=familiesDir+"/"+rank+".faa";
			final ByteFile bf=ByteFile.makeByteFile(file, false);
			String id=null;
			final ByteBuilder seq=new ByteBuilder();
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				if(line.length==0){continue;}
				if(line[0]=='>'){
					if(id!=null){accumulateMember(id, seq.toBytes(), config, famCounts[rank]); famMemberCount[rank]++;}
					id=ProteinSearch.parseHeaderId(line); seq.clear();
				}else{
					seq.append(line);
				}
			}
			if(id!=null){accumulateMember(id, seq.toBytes(), config, famCounts[rank]); famMemberCount[rank]++;}
			bf.close();
			if(famMemberCount[rank]<1){throw new RuntimeException("Family rank "+rank+" ("+file+") has no members.");}
		}
	}

	static void accumulateMember(final String id, final byte[] raw, final AlphaK config, final double[] dst){
		if(ProteinSequence.hasNonterminalStopMarker(raw)){return;}
		final ProteinSequence seq=new ProteinSequence(id, raw);
		addCounts(seq.enc, config.alphabet, config.k, dst);
	}

	/** Runs one (config, profileset, variant) scoring pass over every query against every family. */
	private static Result score(final double[][] qVec, final double[][] famVec, final int[] queryFamily,
			final long[] familyHash, final String[] repIds, final AlphaK config, final String profileSet,
			final String variant, final int[] topNs, final ArrayList<ProteinSequence> queries,
			final int[] famMemberCount, final ByteStreamWriter pqw){
		final int nFamilies=famVec.length;
		final Result out=new Result(config.alphabetName, config.k, config.dims, profileSet, variant, qVec.length, nFamilies, topNs.length);
		final double[] cos=new double[nFamilies];
		final ByteBuilder pb=(pqw==null ? null : new ByteBuilder(256));
		for(int i=0; i<qVec.length; i++){
			final double[] qv=qVec[i];
			double sum=0;
			for(int r=0; r<nFamilies; r++){
				double dot=0;
				final double[] fv=famVec[r];
				for(int d=0; d<qv.length; d++){dot+=qv[d]*fv[d];}
				cos[r]=dot; sum+=dot;
			}
			final int truthFamily=queryFamily[i];
			final double trueCos=cos[truthFamily];
			int best=0, bestOther=-1;
			double bestOtherCos=-Double.MAX_VALUE, maxOtherCos=-Double.MAX_VALUE;
			for(int r=0; r<nFamilies; r++){
				if(r!=truthFamily && cos[r]>maxOtherCos){maxOtherCos=cos[r];}
				if(r!=truthFamily && cos[r]>bestOtherCos){bestOtherCos=cos[r]; bestOther=r;}
				if(betterFamily(r, best, cos, familyHash, repIds)){best=r;}
			}
			int rank=1;
			for(int r=0; r<nFamilies; r++){
				if(r!=truthFamily && betterFamily(r, truthFamily, cos, familyHash, repIds)){rank++;}
			}
			out.trueCosSum+=trueCos;
			out.otherCosSum+=(sum-trueCos)/(nFamilies-1);
			if(best==truthFamily){out.top1++;}
			for(int n=0; n<topNs.length; n++){if(rank<=topNs[n]){out.recall[n]++;}}
			final int famSize=famMemberCount[truthFamily];
			final int qLen=queries.get(i).length();
			out.addStratum(famSizeBucket(famSize), qLenBucket(qLen), rank, topNs);
			if(pb!=null){
				pb.append(config.alphabetName).append('\t').append(config.k).append('\t').append(profileSet).append('\t')
					.append(variant).append('\t').append(queries.get(i).id).append('\t').append(repIds[truthFamily])
					.append('\t').append(rank).append('\t').append(trueCos, 8).append('\t').append(repIds[best])
					.append('\t').append(cos[best], 8).append('\t').append((sum-trueCos)/(nFamilies-1), 8)
					.append('\t').append(maxOtherCos, 8).append('\t').append(qLen).append('\t').append(famSize).nl();
				pqw.print(pb); pb.clear();
			}
		}
		return out;
	}

	private static boolean betterFamily(final int candidate, final int current, final double[] cos,
			final long[] familyHash, final String[] repIds){
		if(cos[candidate]!=cos[current]){return cos[candidate]>cos[current];}
		final int hashCompare=Long.compareUnsigned(familyHash[candidate], familyHash[current]);
		if(hashCompare!=0){return hashCompare<0;}
		return repIds[candidate].compareTo(repIds[current])<0;
	}

	static int famSizeBucket(final int size){return size<100 ? 0 : (size<1000 ? 1 : 2);}
	static int qLenBucket(final int len){return len<100 ? 0 : (len<300 ? 1 : 2);}
	private static final String[] FAMSIZE_LABELS={"<100", "100-999", ">=1000"};
	private static final String[] QLEN_LABELS={"<100", "100-299", ">=300"};

	private static void appendStrata(final ByteStreamWriter stw, final ByteBuilder bb, final Result result,
			final AlphaK config, final String profileSet, final String variant, final int[] topNs){
		for(int b=0; b<3; b++){
			if(result.famSizeN[b]>0){
				bb.append(config.alphabetName).append('\t').append(config.k).append('\t').append(profileSet)
					.append('\t').append(variant).append('\t').append("famsize").append('\t')
					.append(FAMSIZE_LABELS[b]).append('\t').append(result.famSizeN[b]);
				for(int n=0; n<topNs.length; n++){
					bb.append('\t').append(result.famSizeRecall[b][n]/(double)result.famSizeN[b], 8);
				}
				bb.nl(); stw.print(bb); bb.clear();
			}
		}
		for(int b=0; b<3; b++){
			if(result.qLenN[b]>0){
				bb.append(config.alphabetName).append('\t').append(config.k).append('\t').append(profileSet)
					.append('\t').append(variant).append('\t').append("qlen").append('\t')
					.append(QLEN_LABELS[b]).append('\t').append(result.qLenN[b]);
				for(int n=0; n<topNs.length; n++){
					bb.append('\t').append(result.qLenRecall[b][n]/(double)result.qLenN[b], 8);
				}
				bb.nl(); stw.print(bb); bb.clear();
			}
		}
	}

	static String[] loadRankToRepId(final String file, final int maxFamilies){
		final ByteFile bf=ByteFile.makeByteFile(file, true);
		final ArrayList<String> list=new ArrayList<String>();
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0 || line[0]=='#'){continue;}
			final String s=new String(line);
			final int tab=s.indexOf('\t');
			if(tab<0){throw new RuntimeException("Malformed rank_to_repid row: "+s);}
			final int rank=Integer.parseInt(s.substring(0, tab).trim());
			final String rep=s.substring(tab+1).trim();
			if(rank!=list.size()){throw new RuntimeException("rank_to_repid rows must be in rank order starting at 0: expected "+list.size()+", got "+rank);}
			list.add(rep);
			if(maxFamilies>0 && list.size()>=maxFamilies){break;}
		}
		bf.close();
		if(list.isEmpty()){throw new RuntimeException("No rows in rank_to_repid: "+file);}
		return list.toArray(new String[0]);
	}

	private void appendSummaryHeader(final ByteBuilder bb){
		bb.append("alphabet\tk\tdims\tprofileset\tvariant\tqueries\tfamilies\ttop1")
			.append("\tmean_true_cos\tmean_other_cos");
		for(int n : topNs){bb.append("\trecall_at_").append(n);}
		bb.nl();
	}

	private void appendPerQueryHeader(final ByteBuilder bb){
		bb.append("alphabet\tk\tprofileset\tvariant\tquery\ttrue_rep\ttrue_rank\ttrue_cos")
			.append("\ttop1_rep\ttop1_cos\tmean_other_cos\tmax_other_cos\tquery_len\tfamily_size").nl();
	}

	private void appendStrataHeader(final ByteBuilder bb){
		bb.append("alphabet\tk\tprofileset\tvariant\tstratum_type\tstratum_value\tn");
		for(int n : topNs){bb.append("\trecall_at_").append(n);}
		bb.nl();
	}

	static final class Result {
		Result(final String alphabet_, final int k_, final int dims_, final String profileSet_, final String variant_,
				final int queries_, final int families_, final int nCount){
			alphabet=alphabet_; k=k_; dims=dims_; profileSet=profileSet_; variant=variant_;
			queries=queries_; families=families_; recall=new long[nCount];
			famSizeN=new long[3]; qLenN=new long[3];
			famSizeRecall=new long[3][nCount]; qLenRecall=new long[3][nCount];
		}
		void addStratum(final int famSizeB, final int qLenB, final int rank, final int[] topNs){
			famSizeN[famSizeB]++; qLenN[qLenB]++;
			for(int n=0; n<topNs.length; n++){
				if(rank<=topNs[n]){famSizeRecall[famSizeB][n]++; qLenRecall[qLenB][n]++;}
			}
		}
		void append(final ByteBuilder bb, final int[] topNs){
			bb.append(alphabet).append('\t').append(k).append('\t').append(dims).append('\t').append(profileSet)
				.append('\t').append(variant).append('\t').append(queries).append('\t').append(families).append('\t')
				.append(top1/(double)queries, 8).append('\t').append(trueCosSum/queries, 8).append('\t')
				.append(otherCosSum/queries, 8);
			for(long r : recall){bb.append('\t').append(r/(double)queries, 8);}
			bb.nl();
		}
		final String alphabet, profileSet, variant;
		final int k, dims, queries, families;
		long top1;
		double trueCosSum, otherCosSum;
		final long[] recall;
		final long[] famSizeN, qLenN;
		final long[][] famSizeRecall, qLenRecall;
	}

	private static String[] split(final String text){
		final String[] out=text.split(",");
		for(int i=0; i<out.length; i++){out[i]=out[i].trim();}
		return out;
	}
	private static int[] parseInts(final String text){
		final String[] s=split(text); final int[] out=new int[s.length];
		for(int i=0; i<s.length; i++){out[i]=Integer.parseInt(s[i]);}
		return out;
	}

	private String familiesDir, rankToRepIdFile, queryFile, truthFile, outFile, perQueryFile, strataFile;
	private String[] alphabetNames={"aa20", "c14"};
	private int[] ks={1, 2};
	private String[] profileSetNames={"heldout", "insample"};
	private String[] variantNames={"raw", "centered", "logodds"};
	private int[] topNs={1, 5, 25, 50, 100, 1000};
	private double rawPseudocount=0.0;
	private double logPseudocount=1e-3;
	private String backgroundMode="members";
	private int maxFamilies=-1;
}
