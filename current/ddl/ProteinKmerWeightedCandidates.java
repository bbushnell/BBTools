package ddl;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;

import fileIO.ByteFile;
import parse.LineParser1;

/** Exact population-normalized ranking, optionally after complete alphabet aggregation.
 * Used by ProteinKmerCandidates arms=; counts are streamed once per population.
 * @author Yaoyao
 */
public final class ProteinKmerWeightedCandidates{

	static void generate(ProteinKmerCandidates owner, Path real, Path shadow, Path out, String specification,
		int k, String alphabet, long negativeWeight, long positiveWeight, boolean versionedBlacklistHeaders,
		String blacklistProvenance, ArrayList<ProteinKmerCandidates.StableInput> inputs)throws Exception{
		final Arm[] arms=parseArms(specification,k);
		final int[] limits=new int[6];
		for(Arm a:arms){final int at=(a.reduced ? 3 : 0)+(a.weighted ? 1 : a.real ? 2 : 0);limits[at]=Math.max(limits[at],a.size);}
		final int reduced=Math.max(limits[3],Math.max(limits[4],limits[5]));
		require(reduced==0 || alphabet!=null,"Reduced arms require an explicit alphabet");
		final Projection projection=alphabet==null ? null : new Projection(k,alphabet);
		require(projection==null || reduced<=projection.universe,"Requested list exceeds reduced word universe");
		// Score provenance is compressed directly, including small lists, with one consistent backend.
		ProteinKmerCandidates.requireBgzip9();
		final CountSource positive=accepted(owner,real,"real",k,inputs),negative=accepted(owner,shadow,"shadow",k,inputs);
		final Rankings ranks=rankTables(positive,negative,limits,projection,negativeWeight,positiveWeight);
		owner.verifyInputs();Files.createDirectory(out);
		final StringBuilder manifest=new StringBuilder(HEADER).append('\n');
		for(Arm a:arms){
			final Ranking ranked=ranks.get(a.reduced,a.selection);
			final int[] words=Arrays.copyOf(ranked.words,a.size);
			if(a.reduced){for(int i=0;i<words.length;i++){words[i]=projection.representative(words[i]);}}
			Arrays.sort(words);
			final String name=a.name()+".faa"+(ProteinKmerCandidates.fastaBytes(a.size,k)>204800 ? ".gz" : "");
			ProteinKmerCandidates.writeFasta(out.resolve(name),words,k,
				versionedBlacklistHeaders ? blacklistHeader(a, projection, blacklistProvenance) : null);
			manifest.append(name).append('\t').append(k).append('\t').append(a.size).append('\t').append(a.selection)
				.append('\t').append(a.reduced ? "reduced" : "raw").append('\t').append(a.reduced ? projection.alphabet : "none")
				.append('\t').append(a.weighted ? negativeWeight : a.real ? 0 : 1).append('\t').append(a.weighted ? positiveWeight : a.real ? 1 : 0)
				.append('\t').append(negative.windows).append('\t').append(positive.windows).append("\tunmeasured\n");
			final StringBuilder scores=new StringBuilder("rank\tpacked_word\tnegative_occurrences\tpositive_occurrences\tinteger_score\n");
			for(int i=0;i<a.size;i++){
				scores.append(i+1).append('\t').append(a.reduced ? projection.packed(ranked.words[i]) : ranked.words[i])
					.append('\t').append(ranked.negative[i]).append('\t').append(ranked.positive[i]).append('\t').append(ranked.scores[i]).append('\n');
			}
			ProteinKmerCandidates.writeTextArtifact(out.resolve(a.name()+".scores.tsv.gz"),scores.toString());
		}
		final Weights weights=ranks.weights;
		owner.writeText("lists.tsv",manifest.toString());
		owner.writeText("recipe.txt","format=normalized_candidates_v2\n"+(k==6 ? "" : "k="+k+"\n")+"real="+real+"\nshadow="+shadow+"\narms="+specification+
			"\nalphabet="+(projection==null ? "none" : projection.alphabet)+"\nnegative_weight="+negativeWeight+"\npositive_weight="+positiveWeight+
			"\nnegative_valid_windows="+negative.windows+"\npositive_valid_windows="+positive.windows+
			"\nnegative_integer_coefficient="+weights.negative+"\npositive_integer_coefficient="+weights.positive+
			"\nnormalization=each_population_total_valid_windows\naggregation=all_raw_counts_before_reduced_ranking"+
			"\nranking=score_descending_then_packed_word_ascending\nshadow_score=negative_occurrences"+
			"\nweighted_score=negative_occurrences*negative_integer_coefficient+positive_occurrences*positive_integer_coefficient"+
			"\nreduced_fasta="+(projection!=null && projection.dimer ? "first_pair_in_each_ordered_class" : "first_residue_in_each_ordered_class")+"\nraw_rows_real="+positive.rows+"\nraw_rows_shadow="+negative.rows+
			"\nscore_tables=arm_name.scores.tsv.gz\nscore_table_compression=BGZF9"+
			"\nfull_count_top2048_parity=exact_both_populations\nverification=structural\nfile_hashes_computed=0\nadopted=false\n");
		owner.verifyInputs();
		final String pass="PROTEIN_KMER_WEIGHTED_CANDIDATES_PASS lists="+arms.length+" adopted=false";
		owner.writeText("PASS",pass+"\n");System.err.println(pass);
	}

	/** Explicit arms avoid accidental cross-products of sizes and weighting policies. */
	public static Arm[] parseArms(String text){return parseArms(text,6);}
	public static Arm[] parseArms(String text,int k){
		require(k==5 || k==6,"Explicit arms require K5 or K6");
		require(text!=null && !text.isEmpty(),"Require explicit candidate arms");
		final String[] tokens=text.split(",",-1);final Arm[] result=new Arm[tokens.length];final HashSet<String> seen=new HashSet<String>();
		for(int i=0;i<tokens.length;i++){
			final String[] f=tokens[i].split(":",-1);
			require(f.length==3 && (f[0].equals("raw") || f[0].equals("reduced")) &&
				(f[1].equals("shadow") || f[1].equals("weighted") || f[1].equals("real")) && f[2].matches("[0-9]+"),"Invalid arm: "+tokens[i]);
			final int size=Integer.parseInt(f[2]);require(size>0 && size<=ProteinKmerFrequency.universe(k),"Invalid arm size: "+size);
			result[i]=new Arm(f[0].equals("reduced"),f[1],size,k);
			require(seen.add(result[i].name()),"Duplicate arm: "+tokens[i]);
		}return result;
	}
	public static final class Arm{
		private Arm(boolean r,String s,int n,int k_){reduced=r;selection=s;weighted=s.equals("weighted");real=s.equals("real");size=n;k=k_;}
		public String name(){return (reduced ? "reduced" : "raw")+"_"+selection+"_k"+k+"_"+size;}
		public final boolean reduced,weighted,real;public final int size,k;public final String selection;
	}

	private static String blacklistHeader(Arm arm, Projection projection, String provenance){
		require(provenance!=null && safeHeaderValue(provenance),
			"Blacklist provenance must fit one FASTA ID field");
		return "bbgap_blacklist_v1;k="+arm.k+";words="+arm.size+";mode="+
			(arm.reduced ? "reduced_literal" : "literal")+";provenance="+headerValue(provenance)+";aabits="+
			(arm.reduced ? projection.alphabet : "none");
	}

	private static boolean safeHeaderValue(String value){
		for(int i=0; i<value.length(); i++){if(Character.isWhitespace(value.charAt(i))){return false;}}
		return true;
	}

	private static String headerValue(String value){
		final StringBuilder encoded=new StringBuilder(value.length());
		for(int i=0; i<value.length(); i++){
			final char c=value.charAt(i);
			if(c=='%' || c==';'){
				encoded.append('%').append(Character.forDigit((c>>>4)&15, 16)).append(Character.forDigit(c&15, 16));
			}else{encoded.append(c);}
		}
		return encoded.toString();
	}

	private static CountSource accepted(ProteinKmerCandidates owner,Path root,String population,int k,
		ArrayList<ProteinKmerCandidates.StableInput> inputs)throws Exception{
		final Path p=root.resolve(population+"_k"+k),result=p.resolve("result");
		require(!owner.text(p.resolve("accepted.utc")).trim().isEmpty(),"Missing count acceptance");
		for(String name:new String[]{"native.exit","exit","watch.exit"}){require(owner.text(p.resolve(name)).trim().equals("0"),"Count phase failed: "+p);}
		require(owner.text(p.resolve("PASS")).equals("PROTEIN_KMER_COUNTS_PASS\n"),"Missing accepted count PASS");
		require(owner.text(p.resolve("native.time.tsv")).equals(owner.text(p.resolve("accepted.time.tsv"))),"Changed count timing receipt");
		final String sizes=owner.text(p.resolve("accepted.output_sizes.tsv"));
		require(sizes.equals(owner.text(p.resolve("output_sizes.tsv"))),"Changed accepted count sizes");
		owner.bindAcceptedSizes(result,sizes,true);
		require(owner.text(result.resolve("PASS")).equals("PROTEIN_KMER_FREQUENCY_PASS\n"),"Missing count result PASS");
		final String recipe=owner.text(result.resolve("recipe.txt"));
		require(recipe.contains("\npopulation="+population+"\n") && recipe.contains("\nk="+k+"\n") && recipe.contains("\ntop=2048\n"),"Count recipe binding differs");
		final String[] rows=owner.text(result.resolve("summary.tsv")).split("\n",-1);
		require(rows.length==3 && rows[2].isEmpty() && rows[0].equals("population\tk\trecords\tresidues\tvalid_words\tdistinct_words\tselected_words"),"Invalid count summary schema");
		final String[] f=rows[1].split("\t",-1);
		require(f.length==7 && f[0].equals(population) && f[1].equals(Integer.toString(k)) && f[6].equals("2048") &&
			Long.parseLong(f[2])>0 && Long.parseLong(f[3])>0,"Count summary population binding differs");
		final ProteinKmerCandidates.Ranked top=ProteinKmerCandidates.readRankedData(result.resolve("top.tsv.gz"),result.resolve("top.faa.gz"),k,2048,inputs);
		return new CountSource(k,result.resolve("counts.tsv.gz"),Long.parseLong(f[5]),Long.parseLong(f[4]),top);
	}

	/** Fixture/public algorithm entry: sorted full tables, explicit independent population totals. */
	public static Rankings rankTables(Path positive,long positiveRows,long positiveWindows,Path negative,long negativeRows,
		long negativeWindows,int rawLimit,int reducedLimit,String alphabet,long negativeWeight,long positiveWeight)throws Exception{
		return rankTables(new CountSource(6,positive,positiveRows,positiveWindows,null),new CountSource(6,negative,negativeRows,negativeWindows,null),
			new int[]{rawLimit,rawLimit,0,reducedLimit,reducedLimit,0},alphabet==null ? null : new Projection(alphabet),negativeWeight,positiveWeight);
	}
	/** K5/K6 entry includes real-only rankings alongside the historical shadow/weighted choices. */
	public static Rankings rankTables(Path positive,long positiveRows,long positiveWindows,Path negative,long negativeRows,
		long negativeWindows,int rawLimit,int reducedLimit,int k,String alphabet,long negativeWeight,long positiveWeight)throws Exception{
		return rankTables(new CountSource(k,positive,positiveRows,positiveWindows,null),new CountSource(k,negative,negativeRows,negativeWindows,null),
			new int[]{rawLimit,rawLimit,rawLimit,reducedLimit,reducedLimit,reducedLimit},alphabet==null ? null : new Projection(k,alphabet),negativeWeight,positiveWeight);
	}
	private static Rankings rankTables(CountSource positive,CountSource negative,int[] limits,
		Projection projection,long negativeWeight,long positiveWeight)throws Exception{
		final int rawLimit=Math.max(limits[0],Math.max(limits[1],limits[2])),reducedLimit=Math.max(limits[3],Math.max(limits[4],limits[5]));
		for(int limit:limits){require(limit>=0 && limit<=ProteinKmerFrequency.universe(positive.k),"Invalid selection capacity");}
		require(rawLimit>0 || reducedLimit>0,"At least one selection required");
		require(reducedLimit==0 || (projection!=null && reducedLimit<=projection.universe),"Reduced capacity requires sufficient alphabet universe");
		require(positive.k==negative.k && (projection==null || projection.k==positive.k),"Count/projection K differs");
		final Weights weights=new Weights(negative.windows,positive.windows,negativeWeight,positiveWeight);
		final Heap rawShadow=limits[0]==0 ? null : new Heap(limits[0]),rawWeighted=limits[1]==0 ? null : new Heap(limits[1]),rawReal=limits[2]==0 ? null : new Heap(limits[2]);
		final long[] reducedNegative=reducedLimit==0 ? null : new long[projection.universe],
			reducedPositive=reducedLimit==0 ? null : new long[projection.universe];
		try(Cursor p=new Cursor(positive);Cursor n=new Cursor(negative)){
			while(p.valid || n.valid){
				final int word=!n.valid || (p.valid && p.word<n.word) ? p.word : n.word;
				final long pc=p.valid && p.word==word ? p.count : 0,nc=n.valid && n.word==word ? n.count : 0;
				if(rawShadow!=null && nc>0){rawShadow.add(word,nc,nc,pc);}
				if(rawReal!=null && pc>0){rawReal.add(word,pc,nc,pc);}
				if(rawWeighted!=null){rawWeighted.add(word,weights.score(nc,pc),nc,pc);}
				if(reducedLimit>0){final int key=projection.project(word);
					reducedNegative[key]=Math.addExact(reducedNegative[key],nc);reducedPositive[key]=Math.addExact(reducedPositive[key],pc);}
				if(pc>0){p.advance();}if(nc>0){n.advance();}
			}
			p.verify();n.verify();
		}
		Ranking rs=null,rw=null,rr=null;
		if(reducedLimit>0){
			final Heap s=limits[3]==0 ? null : new Heap(limits[3]),w=limits[4]==0 ? null : new Heap(limits[4]),r=limits[5]==0 ? null : new Heap(limits[5]);long ns=0,ps=0;
			for(int key=0;key<projection.universe;key++){
				final long nc=reducedNegative[key],pc=reducedPositive[key];ns=Math.addExact(ns,nc);ps=Math.addExact(ps,pc);
				if(s!=null && nc>0){s.add(key,nc,nc,pc);}if(r!=null && pc>0){r.add(key,pc,nc,pc);}if(w!=null && (nc>0 || pc>0)){w.add(key,weights.score(nc,pc),nc,pc);}
			}
			require(ns==negative.windows && ps==positive.windows,"Reduced aggregation must conserve each population's valid windows");
			rs=s==null ? null : s.finish();rw=w==null ? null : w.finish();rr=r==null ? null : r.finish();
		}
		return new Rankings(rawShadow==null ? null : rawShadow.finish(),rawWeighted==null ? null : rawWeighted.finish(),rawReal==null ? null : rawReal.finish(),rs,rw,rr,weights);
	}

	/** Integer coefficients preserve exact normalized-score ties; overflow fails before reading inputs. */
	public static final class Weights{
		public Weights(long negativeWindows,long positiveWindows,long negativeWeight,long positiveWeight){
			require(negativeWindows>0 && positiveWindows>0 && negativeWeight>0 && positiveWeight>0,"Positive population totals and weights required");
			final long g=gcd(negativeWindows,positiveWindows),wg=gcd(negativeWeight,positiveWeight);
			final long n=Math.multiplyExact(negativeWeight/wg,positiveWindows/g),p=Math.multiplyExact(positiveWeight/wg,negativeWindows/g),cg=gcd(n,p);
			negative=n/cg;positive=p/cg;
			Math.addExact(Math.multiplyExact(negative,negativeWindows),Math.multiplyExact(positive,positiveWindows));
		}
		public long score(long n,long p){
			assert(n>=0 && p>=0) : "Frequency counts cannot be negative; exact weighting assumes valid occurrence counts";
			return Math.addExact(Math.multiplyExact(n,negative),Math.multiplyExact(p,positive));
		}
		public final long negative,positive;
	}
	private static long gcd(long a,long b){assert(a>0 && b>0) : "Weights and denominators must be positive";while(b!=0){final long t=a%b;a=b;b=t;}return a;}

	/** Dense base-alphabet index; sorted index has the same order as packed class symbols. */
	public static final class Projection{
		public Projection(String specification){this(6,specification);}
		public Projection(int k_,String specification){
			require(k_==5 || k_==6,"Projection requires K5 or K6");k=k_;
			require(specification!=null && (specification.startsWith("all:") || specification.startsWith("dimer:")),
				"Require explicit full-word alphabet or nonoverlapping K6 dimer alphabet");
			dimer=specification.startsWith("dimer:");
			require(!dimer || (k==6 && specification.startsWith("dimer:1,3,5:")),
				"Dimer candidate representatives require K6 with nonoverlapping starts1,3,5");
			final ProteinKmerGeometry geometry=ProteinKmerGeometry.make(k,specification);
			if(dimer){
				alphabet="dimer:1,3,5:"+geometry.name().substring("dimer6+1,3,5+".length());
				final String[] groups=alphabet.substring("dimer:1,3,5:".length()).split("/",-1);base=groups.length;
				representatives=new int[base];mapping=null;pairBase=base;bits=32-Integer.numberOfLeadingZeros(base-1);
				universe=Math.multiplyExact(Math.multiplyExact(base,base),base);Arrays.fill(pairs,-1);
				for(int i=0;i<base;i++){
					final String[] members=groups[i].split(",",-1);
					for(int j=0;j<members.length;j++){
						final String pair=members[j];final int raw=(AA.indexOf(pair.charAt(0))<<5)|AA.indexOf(pair.charAt(1));
						pairs[raw]=i;if(j==0){representatives[i]=raw;}
					}
				}
				return;
			}
			alphabet="all:"+geometry.name().substring(("reduced"+k+"+").length());
			final String[] groups=alphabet.substring(4).split("/",-1);base=groups.length;
			representatives=new int[base];mapping=new int[32];Arrays.fill(mapping,-1);
			for(int i=0;i<base;i++){representatives[i]=AA.indexOf(groups[i].charAt(0));for(int j=0;j<groups[i].length();j++){mapping[AA.indexOf(groups[i].charAt(j))]=i;}}
			pairBase=base*base;universe=pairBase*pairBase*(k==6 ? pairBase : base);bits=32-Integer.numberOfLeadingZeros(base-1);
			Arrays.fill(pairs,-1);for(int a=0;a<20;a++){for(int b=0;b<20;b++){pairs[(a<<5)|b]=mapping[a]*base+mapping[b];}}
		}
		public int project(int raw){
			assert(raw>=0 && raw<(1<<(5*k))) : "Raw key must fit K five-bit canonical residues before projection";
			final int a=k==6 ? pairs[raw>>>20] : mapping[raw>>>20],b=pairs[(raw>>>10)&1023],c=pairs[raw&1023];
			assert(a>=0 && b>=0 && c>=0) : "Count reader must reject noncanonical residues before projection";
			return (a*pairBase+b)*pairBase+c;
		}
		public int representative(int key){
			require(key>=0 && key<universe,"Reduced key outside alphabet universe");int raw=0;
			for(int shift=0;shift<5*k;shift+=(dimer ? 10 : 5)){raw|=representatives[key%base]<<shift;key/=base;}return raw;
		}
		public int packed(int key){
			require(key>=0 && key<universe,"Reduced key outside alphabet universe");int packed=0;
			for(int i=0;i<(dimer ? 3 : k);i++){packed|=(key%base)<<(i*bits);key/=base;}return packed;
		}
		public final String alphabet;public final int universe,k;public final boolean dimer;
		private final int base,pairBase,bits;private final int[] representatives,mapping,pairs=new int[1024];
	}

	private static final class CountSource{
		CountSource(int k_,Path p,long r,long w,ProteinKmerCandidates.Ranked t){
			require(k_==5 || k_==6,"Counts require K5 or K6");k=k_;
			require(r>0 && r<=ProteinKmerFrequency.universe(k) && w>=r,"Invalid accepted count totals");path=p;rows=r;windows=w;top=t;
		}
		final int k;final Path path;final long rows,windows;final ProteinKmerCandidates.Ranked top;
	}
	private static final class Cursor implements AutoCloseable{
		Cursor(CountSource source)throws Exception{
			this.source=source;stable=new ProteinKmerCandidates.StableInput(source.path);
			heap=source.top==null ? null : new ProteinKmerCandidates.TopHeap(source.top.words.length);
			reader=ByteFile.makeByteFile(source.path.toString(),false);
			try{
				require(Arrays.equals(reader.nextLine(),"packed_word\tword\toccurrences".getBytes(StandardCharsets.US_ASCII)),"Wrong full-count table schema");advance();
			}catch(Exception | Error e){reader.close();throw e;}
		}
		void advance(){
			final byte[] line=reader.nextLine();valid=line!=null;if(!valid){return;}parser.set(line);
			require(parser.terms()==3 && parser.length(1)==source.k,"Malformed count row for selected K");int key=0;
			for(int i=0;i<source.k;i++){final int b=parser.parseByte(1,i),aa=b>=0 && b<ENCODE.length ? ENCODE[b] : -1;require(aa>=0,"Noncanonical count residue");key=(key<<5)|aa;}
			count=parser.parseLong(2);require(key>word && parser.parseInt(0)==key && count>0,"Invalid count order/key/occurrences");
			word=key;rows++;windows=Math.addExact(windows,count);
			require(rows<=source.rows && windows<=source.windows,"Count table exceeds accepted totals");if(heap!=null){heap.add(key,count);}
		}
		void verify()throws Exception{
			require(!valid && rows==source.rows && windows==source.windows,"Count table does not conserve accepted rows/windows");stable.verify();
			if(heap!=null){final ProteinKmerCandidates.Ranked actual=heap.ranked();
				require(Arrays.equals(actual.words,source.top.words) && Arrays.equals(actual.counts,source.top.counts),"Full counts differ from accepted top2048");}
		}
		@Override public void close()throws Exception{require(!reader.close(),"Count reader failed");stable.verify();}
		boolean valid;int word=-1;long count,rows,windows;final CountSource source;final ByteFile reader;
		final ProteinKmerCandidates.StableInput stable;final ProteinKmerCandidates.TopHeap heap;final LineParser1 parser=new LineParser1('\t');
	}

	public static final class Rankings{
		Rankings(Ranking rs,Ranking rw,Ranking rr,Ranking ds,Ranking dw,Ranking dr,Weights w){rawShadow=rs;rawWeighted=rw;rawReal=rr;reducedShadow=ds;reducedWeighted=dw;reducedReal=dr;weights=w;}
		public Ranking get(boolean reduced,boolean weighted){return reduced ? (weighted ? reducedWeighted : reducedShadow) : (weighted ? rawWeighted : rawShadow);}
		public Ranking get(boolean reduced,String selection){
			require(selection.equals("real") || selection.equals("shadow") || selection.equals("weighted"),"Unknown ranking selection");
			return selection.equals("real") ? (reduced ? reducedReal : rawReal) : get(reduced,selection.equals("weighted"));
		}
		public final Ranking rawShadow,rawWeighted,rawReal,reducedShadow,reducedWeighted,reducedReal;public final Weights weights;
	}
	public static final class Ranking{
		Ranking(int n){words=new int[n];scores=new long[n];negative=new long[n];positive=new long[n];}
		public final int[] words;public final long[] scores,negative,positive;
	}
	private static final class Heap{
		Heap(int capacity){data=new Ranking(capacity);}
		void add(int word,long score,long negative,long positive){
			assert(score>0) : "Only observed words with positive ranking weights enter candidate heap";
			if(size<data.words.length){
				int at=size++;set(at,word,score,negative,positive);
				while(at>0){final int parent=(at-1)/2;if(!worse(at,parent)){break;}swap(at,parent);at=parent;}
			}else if(score>data.scores[0] || (score==data.scores[0] && word<data.words[0])){set(0,word,score,negative,positive);sink();}
		}
		Ranking finish(){
			require(size==data.words.length,"Insufficient distinct observed words for requested blacklist size");final Ranking out=new Ranking(size);
			for(int i=size-1;i>=0;i--){out.words[i]=data.words[0];out.scores[i]=data.scores[0];out.negative[i]=data.negative[0];out.positive[i]=data.positive[0];
				--size;set(0,data.words[size],data.scores[size],data.negative[size],data.positive[size]);if(size>0){sink();}}
			return out;
		}
		private void sink(){int p=0;for(int c=1;c<size;c=p*2+1){if(c+1<size && worse(c+1,c)){c++;}if(!worse(c,p)){break;}swap(c,p);p=c;}}
		private boolean worse(int a,int b){return data.scores[a]<data.scores[b] || (data.scores[a]==data.scores[b] && data.words[a]>data.words[b]);}
		private void set(int at,int w,long s,long n,long p){data.words[at]=w;data.scores[at]=s;data.negative[at]=n;data.positive[at]=p;}
		private void swap(int a,int b){final int w=data.words[a];final long s=data.scores[a],n=data.negative[a],p=data.positive[a];
			set(a,data.words[b],data.scores[b],data.negative[b],data.positive[b]);set(b,w,s,n,p);}
		final Ranking data;int size;
	}
	private static int[] encoding(){final int[] a=new int[128];Arrays.fill(a,-1);for(int i=0;i<AA.length();i++){a[AA.charAt(i)]=i;}return a;}
	private static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	public static final String HEADER="file\tk\tentries\tselection\tspace\talphabet\tnegative_weight\tpositive_weight\tnegative_windows\tpositive_windows\tidentity";
	private static final String AA="ARNDCEQGHILKMFPSTWYV";
	private static final int[] ENCODE=encoding();
}
