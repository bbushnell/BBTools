package prot;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import map.LongHashSet;
import map.LongLongListHashMap;
import parse.LineParser1;
import structures.ByteBuilder;
import structures.IntList;
import structures.LongList;

/**
 * Measures exact reduced-alphabet k-mer shortlist behavior on real protein
 * families. This is the Stage-1B seed-filter assay: it does not score final
 * alignments and must not be interpreted as selecting the authoritative family
 * representation.
 *
 * <p>Each proxy has a proxy id and family id. Shared DISTINCT k-mers are counted
 * per proxy, then multiple proxies are collapsed to a family score by MAX. Equal
 * family scores use unsigned FNV-1a-64 of the UTF-8 family id, then family-id
 * lexical order on the vanishingly unlikely hash collision. Zero-match queries
 * are reported as fallbacks and never credited through arbitrary zero-score tie
 * order.</p>
 *
 * <p>Boundary k-mers keep k symbols: START sentinel + first k-1 residues and
 * last k-1 residues + STOP sentinel. In shared mode START and STOP use the same
 * code; their position within the k-mer disambiguates them, so shared mode
 * requires k&gt;=2. Separate mode uses distinct codes and remains valid at k=1.</p>
 *
 * <p>Usage:
 * {@code java -ea prot.ReducedAlphabetSeedAssay queries=q.faa proxies=p.faa
 * truth=query_to_family.tsv [proxymap=proxy_to_family.tsv]
 * out=seed_assay.tsv [alphabets=aa20,legacy,c6,c7,c8,c9,c12,c14]
 * [ks=4,5,6,7] [boundaries=none,shared,separate] [ns=1,5,10,25]}</p>
 *
 * @author Elly
 */
public final class ReducedAlphabetSeedAssay {

	public static void main(final String[] args){new ReducedAlphabetSeedAssay(args).process();}

	ReducedAlphabetSeedAssay(final String[] args){
		for(String arg : args){
			final int eq=arg.indexOf('=');
			final String a=(eq<0 ? arg : arg.substring(0, eq)).toLowerCase();
			final String b=(eq<0 ? null : arg.substring(eq+1));
			if(a.equals("queries") || a.equals("query")){queryFile=b;}
			else if(a.equals("proxies") || a.equals("refs")){proxyFile=b;}
			else if(a.equals("truth")){truthFile=b;}
			else if(a.equals("proxymap")){proxyMapFile=b;}
			else if(a.equals("out")){outFile=b;}
			else if(a.equals("alphabets")){alphabetNames=split(b);}
			else if(a.equals("ks")){ks=parseInts(b);}
			else if(a.equals("boundaries")){boundaryNames=split(b);}
			else if(a.equals("ns") || a.equals("topn")){topNs=parseInts(b);}
			else if(a.equals("perquery")){perQueryFile=b;}
			else if(a.equals("kmersets")){kmerSetsFile=b;}
			else{throw new RuntimeException("Unknown argument: "+arg);}
		}
		if(queryFile==null || truthFile==null || outFile==null || (proxyFile==null && kmerSetsFile==null)){
			throw new RuntimeException("queries=, truth=, out=, and one of proxies= or kmersets= are required.");
		}
		if(kmerSetsFile!=null && (proxyFile!=null || proxyMapFile!=null)){
			throw new RuntimeException("kmersets= replaces proxies=/proxymap= (a family is scored by its covering set, not by proxy sequences).");
		}
		Arrays.sort(topNs);
		for(int n : topNs){if(n<1){throw new RuntimeException("top-N must be positive: "+n);}}
	}

	private void process(){
		final HashMap<String, String> truth=loadMap(truthFile, "query truth");
		final HashMap<String, String> proxyMap=(proxyMapFile==null ? null :
			loadMap(proxyMapFile, "proxy map"));
		final ArrayList<ProteinSequence> queries=loadFasta(queryFile);
		final Dataset data;
		if(kmerSetsFile==null){
			final ArrayList<ProteinSequence> proxies=loadFasta(proxyFile);
			data=buildDataset(queries, proxies, truth, proxyMap);
		}else{
			//Covering-set mode (plans/COVERING_SET_AA_SPEC_v1.md §4, UMP45 2026-09-02): a family's score is the number of
			//the query's distinct k-mers present in that family's covering set. The file fixes alphabet and k, so the
			//sweep collapses to that one (alphabet,k) with no boundary k-mers; ranking/recall/per-query rows are unchanged.
			final KmerSets sets=loadKmerSets(kmerSetsFile);
			kmerSetAlphabet=sets.alphabet;//may be an explicit symbol list + key, which Alphabet.named() cannot rebuild
			alphabetNames=new String[]{sets.alphabet.name}; ks=new int[]{sets.k}; boundaryNames=new String[]{"none"};
			data=buildKmerSetDataset(queries, sets, truth);
			System.err.println("kmersets: "+sets.families.size()+" families, "+sets.words+" k-mers, alphabet="+
				sets.alphabet.name+" k="+sets.k+" from "+kmerSetsFile);
		}

		final FileFormat ff=FileFormat.testOutput(outFile, FileFormat.TXT, null, true, true, false, false);
		final ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		//Optional per-query distribution output (UMP45, 2026-09-02, Brian's seeding-rate question): one row
		//per (alphabet,k,boundary,query) so the DISTRIBUTION of own-family vs best-other shared k-mers can be
		//plotted, not only the means. Additive: the summary TSV and every existing statistic are unchanged.
		final ByteStreamWriter pqw=(perQueryFile==null ? null :
			new ByteStreamWriter(FileFormat.testOutput(perQueryFile, FileFormat.TXT, null, true, true, false, false)));
		if(pqw!=null){
			pqw.start();
			final ByteBuilder pb=new ByteBuilder(256);
			pb.append("#alphabet\tk\tboundary\tquery\tfamily\tquery_len\tdistinct_kmers\ttrue_score")
				.append("\tbest_other_score\ttrue_rank\tmatched_families\tbest_other_family\tother_score_sum").nl();
			pqw.print(pb);
		}
		try{
			final ByteBuilder bb=new ByteBuilder(4096);
			appendHeader(bb); bsw.print(bb); bb.clear();
			for(String alphabetName : alphabetNames){
				final Alphabet alphabet=(kmerSetAlphabet!=null ? kmerSetAlphabet : Alphabet.named(alphabetName));
				for(int k : ks){
					if(k<1){throw new RuntimeException("k must be positive: "+k);}
					for(String boundaryName : boundaryNames){
						final Boundary boundary=Boundary.named(boundaryName);
						final Result result=run(data, alphabet, k, boundary, topNs, pqw);
						result.append(bb, topNs); bsw.print(bb); bb.clear();
					}
				}
			}
		}finally{
			bsw.poisonAndWait();
			if(pqw!=null){pqw.poisonAndWait();}
		}
	}

	static Result run(final Dataset data, final Alphabet alphabet, final int k,
			final Boundary boundary, final int[] topNs){
		return run(data, alphabet, k, boundary, topNs, null);
	}

	/** As {@link #run(Dataset, Alphabet, int, Boundary, int[])}; when {@code pqw} is non-null, also writes one
	 *  per-query row per scored query (and a row with true_score=0/rank=-1 for zero-match queries) to it. */
	static Result run(final Dataset data, final Alphabet alphabet, final int k,
			final Boundary boundary, final int[] topNs, final ByteStreamWriter pqw){
		final ByteBuilder pb=(pqw==null ? null : new ByteBuilder(256));
		final long buildStart=System.nanoTime();
		final int bits=bitsFor(alphabet.classes+boundary.states);
		if(bits*k>62){throw new RuntimeException("Packed k-mer exceeds 62 bits: alphabet="+
			alphabet.name+", boundary="+boundary.label+", k="+k+", bits="+bits);}
		final LongLongListHashMap index=new LongLongListHashMap();
		long postings=0;
		if(data.proxyWords!=null){
			//Covering-set words were packed at load time with the file's alphabet/k; any other (alphabet,k,boundary)
			//would score them against differently packed query words, so refuse rather than produce numbers.
			if(!alphabet.name.equals(data.proxyAlphabet) || k!=data.proxyK || boundary!=Boundary.NONE){
				throw new RuntimeException("kmersets were packed for alphabet="+data.proxyAlphabet+" k="+data.proxyK+
					" boundary=none; run requested "+alphabet.name+" k="+k+" boundary="+boundary.label);
			}
		}
		for(int proxy=0; proxy<data.proxies.size(); proxy++){
			final LongHashSet kmers=(data.proxyWords!=null ? data.proxyWords[proxy] :
				kmerSet(data.proxies.get(proxy).enc, alphabet, k, boundary));
			for(long word : kmers.toArray()){index.put(word, proxy); postings++;}
		}
		final long buildNanos=System.nanoTime()-buildStart;

		final long queryStart=System.nanoTime();
		final Result out=new Result(alphabet.name, alphabet.classes, k, boundary.label,
			bits, data.queries.size(), data.familyIds.length, data.proxies.size(), topNs.length);
		out.indexKeys=index.size(); out.postings=postings; out.buildNanos=buildNanos;
		final int[] proxyHits=new int[data.proxies.size()];
		final int[] familyScores=new int[data.familyIds.length];
		final IntList touchedProxies=new IntList();
		final IntList touchedFamilies=new IntList();
		final int[] familyQueries=new int[data.familyIds.length];
		final int[] familyScored=new int[data.familyIds.length];
		final int[] familyTop1=new int[data.familyIds.length];
		final long[] familyMarginSum=new long[data.familyIds.length];
		final double[] familyMarginFractionSum=new double[data.familyIds.length];
		final int[][] familyRecall=new int[topNs.length][data.familyIds.length];

		for(int query=0; query<data.queries.size(); query++){
			final int truthFamily=data.queryFamily[query];
			familyQueries[truthFamily]++;
			for(int i=0; i<touchedProxies.size; i++){proxyHits[touchedProxies.get(i)]=0;}
			for(int i=0; i<touchedFamilies.size; i++){familyScores[touchedFamilies.get(i)]=0;}
			touchedProxies.clear(); touchedFamilies.clear();
			final LongHashSet qkmers=kmerSet(data.queries.get(query).enc, alphabet, k, boundary);
			if(qkmers.isEmpty()){out.noValidKmers++; continue;}
			for(long word : qkmers.toArray()){
				final LongList list=index.get(word);
				if(list==null){continue;}
				for(int j=0; j<list.size; j++){
					final int proxy=(int)list.get(j);
					if(proxyHits[proxy]++==0){touchedProxies.add(proxy);}
				}
			}
			for(int i=0; i<touchedProxies.size; i++){
				final int proxy=touchedProxies.get(i), family=data.proxyFamily[proxy];
				if(familyScores[family]==0){touchedFamilies.add(family);}
				if(proxyHits[proxy]>familyScores[family]){familyScores[family]=proxyHits[proxy];}
			}
			if(touchedFamilies.isEmpty()){
				out.zeroMatches++;
				if(pb!=null){
					appendPerQuery(pb, alphabet, k, boundary, data, query, qkmers.size(), 0, 0, -1, 0, -1, 0);
					pqw.print(pb); pb.clear();
				}
				continue;
			}
			out.scoredQueries++;
			familyScored[truthFamily]++;
			out.matchedFamilySum+=touchedFamilies.size;
			int best=-1, bestOtherScore=0, bestOther=-1;
			long otherScoreSum=0;//shared distinct k-mers summed over every family except the truth family (background rate)
			for(int family=0; family<data.familyIds.length; family++){
				if(family!=truthFamily && familyScores[family]>bestOtherScore){
					bestOtherScore=familyScores[family]; bestOther=family;
				}
				if(family!=truthFamily){otherScoreSum+=familyScores[family];}
				if(best<0 || betterFamily(family, best, familyScores, data)){best=family;}
			}
			final int trueScore=familyScores[truthFamily];
			final int margin=trueScore-bestOtherScore;
			out.seedMarginSum+=margin; familyMarginSum[truthFamily]+=margin;
			final double marginFraction=margin/(double)qkmers.size();
			out.seedMarginFractionSum+=marginFraction;
			familyMarginFractionSum[truthFamily]+=marginFraction;
			if(best==truthFamily){out.top1++; familyTop1[truthFamily]++;}
			int rank=-1;
			if(trueScore>0){
				rank=1;
				for(int family=0; family<data.familyIds.length; family++){
					if(family!=truthFamily && betterFamily(family, truthFamily, familyScores, data)){rank++;}
				}
				out.trueRankSum+=rank; out.rankedTruth++;
				for(int n=0; n<topNs.length; n++){
					if(rank<=topNs[n]){out.recall[n]++; familyRecall[n][truthFamily]++;}
				}
			}
			if(pb!=null){
				appendPerQuery(pb, alphabet, k, boundary, data, query, qkmers.size(), trueScore,
					bestOtherScore, rank, touchedFamilies.size, bestOther, otherScoreSum);
				pqw.print(pb); pb.clear();
			}
		}
		out.queryNanos=System.nanoTime()-queryStart;
		int representedFamilies=0, scoredFamilies=0;
		for(int family=0; family<data.familyIds.length; family++){
			if(familyQueries[family]<1){continue;}
			representedFamilies++;
			out.macroTop1+=familyTop1[family]/(double)familyQueries[family];
			if(familyScored[family]>0){
				scoredFamilies++;
				out.macroSeedMargin+=familyMarginSum[family]/(double)familyScored[family];
				out.macroSeedMarginFraction+=familyMarginFractionSum[family]/familyScored[family];
			}
			for(int n=0; n<topNs.length; n++){
				out.macroRecall[n]+=familyRecall[n][family]/(double)familyQueries[family];
			}
		}
		if(representedFamilies>0){
			out.macroTop1/=representedFamilies;
			for(int n=0; n<topNs.length; n++){out.macroRecall[n]/=representedFamilies;}
		}
		if(scoredFamilies>0){
			out.macroSeedMargin/=scoredFamilies;
			out.macroSeedMarginFraction/=scoredFamilies;
		}
		out.queryFamilies=representedFamilies;
		out.scoredFamilies=scoredFamilies;
		return out;
	}

	/** One per-query distribution row; bestOther=-1 / rank=-1 encode "none" (zero-match or unranked truth).
	 *  otherScoreSum = shared distinct k-mers summed over all families except the truth family, so
	 *  otherScoreSum/((families-1)*distinctKmers) is the mean background match rate (Brian's "versus non-consensus"). */
	private static void appendPerQuery(final ByteBuilder pb, final Alphabet alphabet, final int k,
			final Boundary boundary, final Dataset data, final int query, final int distinctKmers,
			final int trueScore, final int bestOtherScore, final int rank, final int matchedFamilies,
			final int bestOther, final long otherScoreSum){
		assert(trueScore>=0 && bestOtherScore>=0 && distinctKmers>=0) : "negative seed count: true="+trueScore+
			" other="+bestOtherScore+" kmers="+distinctKmers+" (scores are counts of shared distinct k-mers, run())";
		assert(trueScore<=distinctKmers && bestOtherScore<=distinctKmers) : "shared k-mers exceed the query's "+
			"distinct k-mers: true="+trueScore+" other="+bestOtherScore+" kmers="+distinctKmers+" (proxyHits counts each "+
			"query k-mer at most once per proxy, run())";
		assert(otherScoreSum>=bestOtherScore && otherScoreSum<=(long)(data.familyIds.length-1)*distinctKmers) :
			"other_score_sum out of range: sum="+otherScoreSum+" best="+bestOtherScore+" families="+data.familyIds.length+
			" kmers="+distinctKmers+" (sum over the other families of a per-family MAX over proxies, each <= distinctKmers, run())";
		final ProteinSequence q=data.queries.get(query);
		pb.append(alphabet.name).append('\t').append(k).append('\t').append(boundary.label).append('\t')
			.append(q.id).append('\t').append(data.familyIds[data.queryFamily[query]]).append('\t')
			.append(q.length()).append('\t').append(distinctKmers).append('\t').append(trueScore).append('\t')
			.append(bestOtherScore).append('\t').append(rank).append('\t').append(matchedFamilies).append('\t')
			.append(bestOther<0 ? "-" : data.familyIds[bestOther]).append('\t').append(otherScoreSum).nl();
	}

	private static boolean betterFamily(final int candidate, final int current,
			final int[] score, final Dataset data){
		if(score[candidate]!=score[current]){return score[candidate]>score[current];}
		final int hashCompare=Long.compareUnsigned(data.familyHash[candidate], data.familyHash[current]);
		if(hashCompare!=0){return hashCompare<0;}
		return data.familyIds[candidate].compareTo(data.familyIds[current])<0;
	}

	static LongHashSet kmerSet(final byte[] enc, final Alphabet alphabet, final int k,
			final Boundary boundary){
		final LongHashSet set=new LongHashSet(Math.max(2, enc.length+2));
		fillKmerSetInternal(enc,alphabet,k,boundary,set,null);
		return set;
	}

	/**
	 * Fills caller-owned scratch with every distinct packed k-mer in encounter order.  This is the
	 * allocation-free sibling of {@link #kmerSet}: both entry points share the same implementation,
	 * so shortlist scoring and direct-family calibration cannot drift in encoding, ambiguous-residue
	 * reset, boundary handling, or duplicate suppression.  Both destinations are cleared first.
	 *
	 * @param enc Encoded protein residues.
	 * @param alphabet Reduced alphabet used to pack residues.
	 * @param k K-mer length.
	 * @param boundary Boundary-sentinel policy.
	 * @param set Reused membership scratch.
	 * @param words Reused encounter-order list of the values newly inserted into {@code set}.
	 */
	static void fillKmerSet(final byte[] enc, final Alphabet alphabet, final int k,
			final Boundary boundary, final LongHashSet set, final LongList words){
		if(set==null || words==null){throw new IllegalArgumentException("fillKmerSet requires non-null caller-owned set and word list.");}
		fillKmerSetInternal(enc,alphabet,k,boundary,set,words);
	}

	/** Shared implementation; {@code words} is null only for the allocating legacy entry point. */
	private static void fillKmerSetInternal(final byte[] enc, final Alphabet alphabet, final int k,
			final Boundary boundary, final LongHashSet set, final LongList words){
		if(enc==null || alphabet==null || boundary==null || set==null || k<1){
			throw new IllegalArgumentException("Invalid k-mer fill arguments: enc="+(enc==null ? "null" : enc.length)+
				" alphabet="+(alphabet==null ? "null" : alphabet.name)+" k="+k+" boundary="+boundary+" set="+(set!=null)+".");
		}
		set.clear();
		if(words!=null){words.clear();}
		if(boundary==Boundary.SHARED && k<2){
			throw new RuntimeException("Shared START/STOP boundary mode requires k>=2 so "
				+"sentinel position can disambiguate the two boundaries.");
		}
		final int bits=bitsFor(alphabet.classes+boundary.states);
		if(bits*k>62){throw new RuntimeException("k-mer needs "+(bits*k)+" bits.");}
		final long mask=(1L<<(bits*k))-1;
		long word=0;
		int valid=0;
		for(byte residue : enc){
			final int code=alphabet.code(residue);
			if(code<0){word=0; valid=0; continue;}
			word=((word<<bits)|code)&mask;
			if(++valid>=k && set.add(word) && words!=null){words.add(word);}
		}
		if(boundary.states>0 && enc.length>=k-1){
			long start=boundary.startCode(alphabet.classes);
			boolean validStart=true;
			for(int i=0; i<k-1; i++){
				final int code=alphabet.code(enc[i]);
				if(code<0){validStart=false; break;}
				start=(start<<bits)|code;
			}
			if(validStart && set.add(start) && words!=null){words.add(start);}
			long stop=0;
			boolean validStop=true;
			for(int i=enc.length-(k-1); i<enc.length; i++){
				final int code=alphabet.code(enc[i]);
				if(code<0){validStop=false; break;}
				stop=(stop<<bits)|code;
			}
			final long stopWord=(stop<<bits)|boundary.stopCode(alphabet.classes);
			if(validStop && set.add(stopWord) && words!=null){words.add(stopWord);}
		}
	}

	static int bitsFor(final int states){
		return ReducedAlphabet.bitsFor(states);
	}

	static long stableHash(final String id){
		long hash=0xcbf29ce484222325L;
		for(byte b : id.getBytes(StandardCharsets.UTF_8)){
			hash^=(b&0xFF); hash*=0x100000001b3L;
		}
		return hash;
	}

	static Dataset buildDataset(final ArrayList<ProteinSequence> queries,
			final ArrayList<ProteinSequence> proxies, final HashMap<String, String> truth,
			final HashMap<String, String> proxyMap){
		if(queries.isEmpty()){throw new RuntimeException("Query FASTA contains no non-empty sequences.");}
		if(proxies.isEmpty()){throw new RuntimeException("Proxy FASTA contains no non-empty sequences.");}
		final HashSet<String> ids=new HashSet<String>();
		for(ProteinSequence sequence : queries){
			if(!ids.add(sequence.id)){throw new RuntimeException("Duplicate query id: "+sequence.id);}
			if(!truth.containsKey(sequence.id)){throw new RuntimeException("Missing truth for query: "+sequence.id);}
		}
		ids.clear();
		final HashSet<String> families=new HashSet<String>();
		final String[] proxyFamilyIds=new String[proxies.size()];
		for(int i=0; i<proxies.size(); i++){
			final ProteinSequence proxy=proxies.get(i);
			if(!ids.add(proxy.id)){throw new RuntimeException("Duplicate proxy id: "+proxy.id);}
			final String family=(proxyMap==null ? proxy.id : proxyMap.get(proxy.id));
			if(family==null){throw new RuntimeException("Missing family mapping for proxy: "+proxy.id);}
			proxyFamilyIds[i]=family; families.add(family);
		}
		final String[] familyIds=families.toArray(new String[0]);
		Arrays.sort(familyIds);
		final HashMap<String, Integer> familyIndex=new HashMap<String, Integer>(familyIds.length*2);
		final long[] familyHash=new long[familyIds.length];
		for(int i=0; i<familyIds.length; i++){familyIndex.put(familyIds[i], i); familyHash[i]=stableHash(familyIds[i]);}
		final int[] proxyFamily=new int[proxies.size()];
		for(int i=0; i<proxyFamily.length; i++){proxyFamily[i]=familyIndex.get(proxyFamilyIds[i]);}
		final int[] queryFamily=new int[queries.size()];
		for(int i=0; i<queries.size(); i++){
			final String family=truth.get(queries.get(i).id);
			final Integer index=familyIndex.get(family);
			if(index==null){throw new RuntimeException("Truth family has no proxy: query="+
				queries.get(i).id+", family="+family);}
			queryFamily[i]=index;
		}
		return new Dataset(queries, proxies, queryFamily, proxyFamily, familyIds, familyHash);
	}

	/**
	 * Loads a covering-set table (spec §3.6): comment lines {@code #alphabet=<name>} and {@code #k=<int>} fix the
	 * packing; data rows {@code family_id<TAB>kmer<TAB>...} give one k-mer per row in the alphabet's letters (for a
	 * reduced alphabet, the group-representative letters). Every k-mer must have exactly k letters that encode under
	 * the alphabet; a family listing the same k-mer twice is a contract violation and fails loud.
	 */
	static KmerSets loadKmerSets(final String file){
		final ByteFile bf=ByteFile.makeByteFile(file, true);
		final LineParser1 lp=new LineParser1((byte)'\t');
		//Header contract (spec §3.6 + CoveringSet's writer, Sayu a7aeafa): either "#alphabet=<name>" / "#k=<int>", or the
		//writer's tab form "#alphabet<TAB><symbols>", "#key<TAB><groups|none>", "#k<TAB><k>[<TAB>kdesign ...]". With a key the
		//alphabet is the key's groups (symbols = group representatives); without one, an explicit symbol list is a
		//singleton-group alphabet and a known name (aa20/amino/legacy/c8/c14/...) is that preset.
		String alphabetName=null, keySpec=null; int k=-1; Alphabet alphabet=null; int bits=0; long mask=0;
		final HashMap<String, LongHashSet> families=new HashMap<String, LongHashSet>();
		final ArrayList<String> order=new ArrayList<String>();
		long words=0, row=0;
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			row++;
			if(line.length==0){continue;}
			if(line[0]=='#'){
				final String s=new String(line);
				if(s.startsWith("#alphabet=")){alphabetName=s.substring(10).trim();}
				else if(s.startsWith("#alphabet\t")){alphabetName=s.substring(10).trim();}
				else if(s.startsWith("#key\t")){keySpec=s.substring(5).trim(); if(keySpec.equalsIgnoreCase("none")){keySpec=null;}}
				else if(s.startsWith("#k=")){k=Integer.parseInt(s.substring(3).trim());}
				else if(s.startsWith("#k\t")){k=Integer.parseInt(s.substring(3).split("\t")[0].trim());}
				continue;
			}
			if(alphabet==null){
				if(alphabetName==null || k<1){throw new RuntimeException("kmersets file "+file+" must declare #alphabet and #k before data rows (row "+row+")");}
				alphabet=alphabetFor(alphabetName, keySpec);
				bits=bitsFor(alphabet.classes);
				if(bits*k>62){throw new RuntimeException("kmersets: "+bits+" bits x k="+k+" exceeds 62 bits");}
				mask=(1L<<(bits*k))-1;
			}
			lp.set(line);
			if(lp.terms()<2){throw new RuntimeException("Malformed kmersets row "+row+": "+new String(line));}
			final String family=lp.parseString(0);
			final int len=lp.setBounds(1);//LineParser1.setBounds(term) selects term 1; a()/b() are its [start,end)
			final int a=lp.a(), b=lp.b();
			if(len!=k){throw new RuntimeException("kmersets row "+row+": k-mer length "+len+" != k="+k+" for family "+family);}
			long word=0;
			for(int i=a; i<b; i++){
				final int c=line[i]&0xFF;
				final int encoded=(c<AminoAcid.acidToNumberExtended.length ? AminoAcid.acidToNumberExtended[c] : -1);
				final int code=(encoded<0 ? -1 : alphabet.code((byte)encoded));
				if(code<0){throw new RuntimeException("kmersets row "+row+": residue '"+(char)c+"' is not in alphabet "+alphabet.name);}
				word=((word<<bits)|code)&mask;
			}
			LongHashSet set=families.get(family);
			if(set==null){set=new LongHashSet(64); families.put(family, set); order.add(family);}
			if(!set.add(word)){throw new RuntimeException("kmersets row "+row+": duplicate k-mer in family "+family);}
			words++;
		}
		bf.close();
		if(families.isEmpty()){throw new RuntimeException("kmersets file has no data rows: "+file);}
		return new KmerSets(alphabet, k, families, order, words);
	}

	/** Alphabet for a covering-set file: a key (slash groups) defines it outright; otherwise a preset name, or an explicit
	 *  symbol list as singleton groups (each letter its own class -- only meaningful when the letters are real residues). */
	static Alphabet alphabetFor(final String alphabetSpec, final String keySpec){
		if(keySpec!=null && keySpec.length()>0){
			final String groups=(keySpec.indexOf('/')>=0 ? keySpec.toUpperCase() : null);
			if(groups==null){return Alphabet.named(keySpec);}//a preset name used as the key
			return new Alphabet(alphabetSpec, groups);
		}
		final String name=alphabetSpec.trim();
		if(name.equalsIgnoreCase("amino")){return Alphabet.named("aa20");}
		try{return Alphabet.named(name);}catch(RuntimeException e){/* not a preset: explicit symbol list */}
		final StringBuilder sb=new StringBuilder(name.length()*2);
		for(int i=0; i<name.length(); i++){if(i>0){sb.append('/');} sb.append(Character.toUpperCase(name.charAt(i)));}
		return new Alphabet(name, sb.toString());
	}

	/** Dataset for covering-set mode: one placeholder proxy per family (id = family id) whose word set is the family's
	 *  covering set; the truth map's families must all have a set. Family index order = sorted family ids, as in buildDataset. */
	static Dataset buildKmerSetDataset(final ArrayList<ProteinSequence> queries, final KmerSets sets,
			final HashMap<String, String> truth){
		if(queries.isEmpty()){throw new RuntimeException("Query FASTA contains no non-empty sequences.");}
		final HashSet<String> ids=new HashSet<String>();
		for(ProteinSequence sequence : queries){
			if(!ids.add(sequence.id)){throw new RuntimeException("Duplicate query id: "+sequence.id);}
			if(!truth.containsKey(sequence.id)){throw new RuntimeException("Missing truth for query: "+sequence.id);}
		}
		final String[] familyIds=sets.order.toArray(new String[0]);
		Arrays.sort(familyIds);
		final ArrayList<ProteinSequence> proxies=new ArrayList<ProteinSequence>(familyIds.length);
		final LongHashSet[] proxyWords=new LongHashSet[familyIds.length];
		final int[] proxyFamily=new int[familyIds.length];
		final long[] familyHash=new long[familyIds.length];
		final HashMap<String, Integer> familyIndex=new HashMap<String, Integer>(familyIds.length*2);
		for(int i=0; i<familyIds.length; i++){
			proxies.add(new ProteinSequence(familyIds[i], "A"));//placeholder: only the id is used in this mode
			proxyWords[i]=sets.families.get(familyIds[i]);
			proxyFamily[i]=i; familyHash[i]=stableHash(familyIds[i]); familyIndex.put(familyIds[i], i);
		}
		final int[] queryFamily=new int[queries.size()];
		for(int i=0; i<queries.size(); i++){
			final String family=truth.get(queries.get(i).id);
			final Integer index=familyIndex.get(family);
			if(index==null){throw new RuntimeException("Truth family has no covering set: query="+queries.get(i).id+", family="+family);}
			queryFamily[i]=index;
		}
		return new Dataset(queries, proxies, queryFamily, proxyFamily, familyIds, familyHash,
			proxyWords, sets.alphabet.name, sets.k);
	}

	static final class KmerSets {
		KmerSets(Alphabet alphabet_, int k_, HashMap<String, LongHashSet> families_, ArrayList<String> order_, long words_){
			alphabet=alphabet_; k=k_; families=families_; order=order_; words=words_;
		}
		final Alphabet alphabet; final int k; final HashMap<String, LongHashSet> families; final ArrayList<String> order; final long words;
	}

	static HashMap<String, String> loadMap(final String file, final String label){
		final HashMap<String, String> map=new HashMap<String, String>();
		final ByteFile bf=ByteFile.makeByteFile(file, true);
		final LineParser1 lp=new LineParser1((byte)'\t');
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0 || line[0]=='#'){continue;}
			lp.set(line);
			if(lp.terms()<2){throw new RuntimeException("Malformed "+label+" row: "+new String(line));}
			final String key=lp.parseString(0), value=lp.parseString(1);
			if(map.put(key, value)!=null){throw new RuntimeException("Duplicate "+label+" key: "+key);}
		}
		bf.close();
		return map;
	}

	/**
	 * Loads a protein FASTA through the assigner's own loader ({@link ProteinSearch#readFasta}: ByteFile lines,
	 * first-whitespace-token id, one terminal '*' stripped, nonterminal-'*' records skipped and counted), so the
	 * assay scores exactly the sequences {@link ProteinSearcher} would see. The previous path went through
	 * {@code ConcurrentReadInputStream}, whose FASTA reader takes amino mode from the file EXTENSION
	 * ({@code FileFormat} line 363, {@code amino=isAminoExt(rawExtension())}); a {@code .fasta} proxy file was
	 * therefore validated as nucleotides and rejected at the first 'I' (Dori job 25480758, UMP45 2026-09-02).
	 */
	static ArrayList<ProteinSequence> loadFasta(final String file){
		final ArrayList<ProteinSequence> out=new ArrayList<ProteinSequence>(ProteinSearch.readFasta(file));
		assert(!out.isEmpty()) : "ProteinSearch.readFasta throws on an empty FASTA, so a returned list is non-empty: "+file;
		return out;
	}

	private void appendHeader(final ByteBuilder bb){
		bb.append("alphabet\tclasses\tk\tboundary\tbits\tqueries\tscored_queries")
			.append("\tno_valid_kmers\tzero_matches\tfamilies\tproxies\tindex_keys\tpostings")
			.append("\tquery_families\tscored_families\ttop1_micro\ttop1_macro")
			.append("\tmean_seed_margin_scored_micro\tmean_seed_margin_scored_macro")
			.append("\tmean_seed_margin_fraction_scored_micro\tmean_seed_margin_fraction_scored_macro")
			.append("\tmean_true_rank_matched_truth\tmean_matched_families")
			.append("\tbuild_seconds\tquery_seconds");
		for(int n : topNs){bb.append("\trecall_at_").append(n).append("_micro");}
		for(int n : topNs){bb.append("\trecall_at_").append(n).append("_macro");}
		bb.nl();
	}

	static final class Dataset {
		Dataset(ArrayList<ProteinSequence> queries_, ArrayList<ProteinSequence> proxies_,
				int[] queryFamily_, int[] proxyFamily_, String[] familyIds_, long[] familyHash_){
			this(queries_, proxies_, queryFamily_, proxyFamily_, familyIds_, familyHash_, null, null, -1);
		}
		/** proxyWords non-null = covering-set mode: proxy i's word set is given (packed for proxyAlphabet/proxyK, no boundary). */
		Dataset(ArrayList<ProteinSequence> queries_, ArrayList<ProteinSequence> proxies_,
				int[] queryFamily_, int[] proxyFamily_, String[] familyIds_, long[] familyHash_,
				LongHashSet[] proxyWords_, String proxyAlphabet_, int proxyK_){
			queries=queries_; proxies=proxies_; queryFamily=queryFamily_; proxyFamily=proxyFamily_;
			familyIds=familyIds_; familyHash=familyHash_;
			proxyWords=proxyWords_; proxyAlphabet=proxyAlphabet_; proxyK=proxyK_;
			assert(proxyWords==null || proxyWords.length==proxies.size()) : "one word set per proxy: "+
				(proxyWords==null ? -1 : proxyWords.length)+" vs "+proxies.size()+" (buildKmerSetDataset)";
		}
		final ArrayList<ProteinSequence> queries, proxies;
		final int[] queryFamily, proxyFamily;
		final String[] familyIds;
		final long[] familyHash;
		final LongHashSet[] proxyWords;
		final String proxyAlphabet;
		final int proxyK;
	}

	static final class Alphabet {
		Alphabet(final String name_, final String groups){
			name=name_; delegate=ReducedAlphabet.parse(representatives(groups), groups);
			classes=delegate.classes();
		}
		int code(final byte encoded){return delegate.codeEncoded(encoded);}
		static Alphabet named(final String raw){
			final String name=raw.toLowerCase();
			final ReducedAlphabet shared=ReducedAlphabet.named(name);
			return new Alphabet(shared);
		}
		private static String representatives(final String groups){
			final String[] tokens=groups.split("/", -1); final StringBuilder sb=new StringBuilder(tokens.length);
			for(String token : tokens){if(token.length()==0){throw new RuntimeException("Empty alphabet group");} sb.append(token.charAt(0));}
			return sb.toString();
		}
		private Alphabet(final ReducedAlphabet shared){name=shared.name(); delegate=shared; classes=shared.classes();}
		final String name;
		final int classes;
		final ReducedAlphabet delegate;
	}

	enum Boundary {
		NONE("none", 0), SHARED("shared", 1), SEPARATE("separate", 2);
		Boundary(String label_, int states_){label=label_; states=states_;}
		static Boundary named(final String raw){
			for(Boundary b : values()){if(b.label.equalsIgnoreCase(raw)){return b;}}
			throw new RuntimeException("Unknown boundary mode: "+raw);
		}
		int startCode(final int classes){return classes;}
		int stopCode(final int classes){return this==SEPARATE ? classes+1 : classes;}
		final String label;
		final int states;
	}

	static final class Result {
		Result(String alphabet_, int classes_, int k_, String boundary_, int bits_,
				int queries_, int families_, int proxies_, int nCount){
			alphabet=alphabet_; classes=classes_; k=k_; boundary=boundary_; bits=bits_;
			queries=queries_; families=families_; proxies=proxies_;
			recall=new long[nCount]; macroRecall=new double[nCount];
		}
		void append(final ByteBuilder bb, final int[] topNs){
			bb.append(alphabet).append('\t').append(classes).append('\t').append(k).append('\t')
				.append(boundary).append('\t').append(bits).append('\t').append(queries).append('\t')
				.append(scoredQueries).append('\t').append(noValidKmers).append('\t').append(zeroMatches)
				.append('\t').append(families).append('\t').append(proxies).append('\t').append(indexKeys)
				.append('\t').append(postings).append('\t').append(queryFamilies).append('\t')
				.append(scoredFamilies).append('\t').append(ratio(top1, queries), 8)
				.append('\t').append(macroTop1, 8).append('\t').append(ratio(seedMarginSum, scoredQueries), 8)
				.append('\t').append(macroSeedMargin, 8).append('\t')
				.append(scoredQueries<1 ? Double.NaN : seedMarginFractionSum/scoredQueries, 8)
				.append('\t').append(macroSeedMarginFraction, 8)
				.append('\t').append(ratio(trueRankSum, rankedTruth), 8)
				.append('\t').append(ratio(matchedFamilySum, scoredQueries), 8)
				.append('\t').append(buildNanos/1e9, 6).append('\t').append(queryNanos/1e9, 6);
			for(long count : recall){bb.append('\t').append(ratio(count, queries), 8);}
			for(double value : macroRecall){bb.append('\t').append(value, 8);}
			bb.nl();
		}
		private static double ratio(final long numerator, final long denominator){
			return denominator<1 ? Double.NaN : numerator/(double)denominator;
		}
		final String alphabet, boundary;
		final int classes, k, bits, queries, families, proxies;
		long indexKeys, postings, scoredQueries, noValidKmers, zeroMatches, top1;
		int queryFamilies, scoredFamilies;
		long seedMarginSum, trueRankSum, rankedTruth, matchedFamilySum, buildNanos, queryNanos;
		final long[] recall;
		double seedMarginFractionSum, macroTop1, macroSeedMargin, macroSeedMarginFraction;
		final double[] macroRecall;
	}

	private static String[] split(final String text){
		final String[] out=text.split(",");
		for(int i=0; i<out.length; i++){out[i]=out[i].trim();}
		return out;
	}
	private static int[] parseInts(final String text){
		final String[] split=split(text); final int[] out=new int[split.length];
		for(int i=0; i<split.length; i++){out[i]=Integer.parseInt(split[i]);}
		return out;
	}

	private String queryFile, proxyFile, truthFile, proxyMapFile, outFile, perQueryFile, kmerSetsFile;
	private Alphabet kmerSetAlphabet;
	private String[] alphabetNames={"aa20", "legacy", "c6", "c7", "c8", "c9", "c12", "c14"};
	private int[] ks={4, 5, 6, 7};
	private String[] boundaryNames={"none", "shared", "separate"};
	private int[] topNs={1, 5, 10, 25};
}
