package ddl;

import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.LinkOption;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.nio.file.attribute.BasicFileAttributes;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Locale;
import java.util.Objects;

import dna.Data;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import map.IntHashSet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;

/** Deterministic comparison lists from separately accepted real/shadow frequency counts.
 * This tool creates candidates only; it never installs a filter or measures model accuracy.
 * @author Yaoyao
 */
public final class ProteinKmerCandidates{

	public static void main(String[] args) throws Exception{
		final ProteinKmerCandidates tool=new ProteinKmerCandidates(args);
		try{tool.process();}finally{Shared.closeStream(tool.outstream);}
	}

	public ProteinKmerCandidates(String[] args){
		final PreParser pp=new PreParser(args,getClass(),false);outstream=pp.outstream;
		final Parser parser=new Parser();parser.overwrite=false;
		final HashSet<String> seen=new HashSet<String>();
		for(String arg:pp.args){
			final String[] pair=Parse.splitOnFirst(arg,'=');
			final String key=pair[0].toLowerCase(Locale.ROOT),value=pair.length>1 ? pair[1] : null;
			require(seen.add(key),"Duplicate option: "+key);
			if(key.equals("real")){realRoot=Paths.get(value);}
			else if(key.equals("rankedexport")){rankedExport=Paths.get(value);}
			else if(key.equals("shadow")){shadowRoot=Paths.get(value);}
			else if(key.equals("legacy")){legacy=Paths.get(value);}
			else if(key.equals("arms")){arms=value;}
			else if(key.equals("alphabet")){alphabet=value;}
			else if(key.equals("blacklistprovenance")){blacklistProvenance=value;}
			else if(key.equals("versionedblacklistheaders")){versionedBlacklistHeaders=Parse.parseBoolean(value);}
			else if(key.equals("k")){candidateK=Integer.parseInt(value);}
			else if(key.equals("negativeweight")){negativeWeight=Long.parseLong(value);}
			else if(key.equals("positiveweight")){positiveWeight=Long.parseLong(value);}
			else if(key.equals("k6sizes")){
				k6sizes=parseSizes(value);
			}
			else if(key.equals("t") || key.equals("threads")){require(Integer.parseInt(value)==1,"Use t=1");}
			else if(!parser.parse(arg,key,value)){throw new IllegalArgumentException("Unknown option: "+arg);}
		}
		if(rankedExport!=null){
			require(parser.out1!=null && realRoot==null && shadowRoot==null && legacy==null && arms==null && alphabet==null &&
				!seen.contains("k") && !seen.contains("k6sizes") && !seen.contains("negativeweight") && !seen.contains("positiveweight") &&
				!seen.contains("versionedblacklistheaders") && !seen.contains("blacklistprovenance"),
				"rankedexport requires only an accepted candidate root and fresh out; no selection parameters");
		}else{
			require(realRoot!=null && shadowRoot!=null && (legacy!=null || arms!=null) && parser.out1!=null,"Require real, shadow, out and either legacy or arms");
			require(arms==null ? !seen.contains("k") && alphabet==null && negativeWeight==4 && positiveWeight==1 : legacy==null && !seen.contains("k6sizes"),
				"Explicit arms replace legacy/k6sizes; k, alphabet and weights require arms");
			require(arms!=null || (!seen.contains("versionedblacklistheaders") && !seen.contains("blacklistprovenance")),
				"Versioned blacklist headers require explicit arms");
		}
		require(candidateK==5 || candidateK==6,"Explicit-arm k must be 5 or 6");
		require(parser.in1==null && parser.in2==null && parser.out2==null && !parser.append && !parser.overwrite && parser.maxReads<0,
			"Require complete inputs and fresh output directory");
		output=Paths.get(parser.out1);
		maxK6=k6sizes[k6sizes.length-1];
		for(Path p:new Path[]{realRoot,shadowRoot,legacy,rankedExport,output}){require(p==null || p.isAbsolute(),"Require absolute paths: "+p);}
	}

	private void process() throws Exception{
		require(!Files.exists(output,LinkOption.NOFOLLOW_LINKS),"Output already exists: "+output);
		Shared.setThreads(1);Shared.AMINO_IN=true;ReadWrite.setZipThreads(1);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=false;
		if(rankedExport!=null){exportRanked();return;}
		if(arms!=null){
			ProteinKmerWeightedCandidates.generate(this, realRoot, shadowRoot, output, arms, candidateK, alphabet, negativeWeight,
				positiveWeight, versionedBlacklistHeaders, blacklistProvenance, inputs);
			return;
		}
		final boolean compressed=fastaBytes(maxK6,6)>204800;
		if(compressed){requireBgzip9();}
		final ArrayList<Selection> lists=new ArrayList<Selection>();
		final StringBuilder recipe=new StringBuilder("real="+realRoot+"\nshadow="+shadowRoot+"\nlegacy="+legacy+
			"\ntarget_entries=2048\nk6sizes="+sizesText()+"\nword_semantics=literal_raw_amino_acids\nalphabet="+AA+
			"\nsource_order=occurrences_descending_then_packed_ascending\n"+
			"mixed_selection=real_rank1,shadow_rank1,real_rank2,shadow_rank2,...\n"+
			"duplicate_policy=advance_source_cursor_on_every_encounter\nstop=immediately_on_requested_unique_count\n"+
			"union=full_union_different_size_control\nemission_order=ascending_canonical_five_bit_packed_word\n"+
			"legacy_emission=same_word_set_canonical_order_new_headers\nidentity=unmeasured\nfile_hashes_computed=0\n"+
			"verification=structural_regular_file_canonical_target_size_mtime_filekey\n"+
			"acceptance_binding=consumed_result_sizes_match_accepted_receipt\n"+
			"full_count_selection="+(maxK6>LIMIT ? "bounded_primitive_heap_top"+maxK6 : "not_read")+"\n"+
			"large_fasta_compression=BGZF9_if_raw_bytes_exceed204800\n"+
			"large_fasta_backend="+(compressed ? "external_bgzip_-c_-f_-l_9_-@_1" : "not_used")+"\n"+
			"bgzip_version="+(compressed ? Data.BGZIP_VERSION : "not_used")+"\n"+
			"cross_population_frequency_coverage=not_measured\nadopted=false\n");
		for(int k:new int[]{5,6}){
			final int[] real=accepted(realRoot,"real",k,recipe),shadow=accepted(shadowRoot,"shadow",k,recipe);
			final int[] r=sortedUnique(Arrays.copyOf(real,LIMIT)),s=sortedUnique(Arrays.copyOf(shadow,LIMIT));
			lists.add(new Selection("real",k,r,r,s));lists.add(new Selection("shadow",k,s,r,s));
			lists.add(new Selection("mixed",k,mixed(real,shadow,LIMIT),r,s));
			lists.add(new Selection("union",k,union(r,s),r,s));
			if(k==6){
				int[] previous=new int[0];
				for(int size:k6sizes){
					final int[] selected=mixed(real,shadow,size);
					for(int word:previous){require(Arrays.binarySearch(selected,word)>=0,"Mixed candidates are not nested");}
					previous=selected;
					if(size==LIMIT){continue;}
					lists.add(new Selection("real_"+size,k,sortedUnique(Arrays.copyOf(real,size)),r,s));
					lists.add(new Selection("shadow_"+size,k,sortedUnique(Arrays.copyOf(shadow,size)),r,s));
					lists.add(new Selection("mixed_"+size,k,selected,r,s));
				}
			}
		}
		final int[] old=readFasta(legacy,5,LIMIT,null,null,inputs);
		lists.add(new Selection("legacy",5,sortedUnique(old),lists.get(0).real,lists.get(0).shadow));
		verifyInputs();Files.createDirectory(output);
		final StringBuilder summary=new StringBuilder("file\tk\tentries\tselection\tin_real_top2048\tin_shadow_top2048\tin_both_top2048\tidentity\n");
		for(Selection list:lists){
			final String name=list.name+".k"+list.k+".faa"+(fastaBytes(list.words.length,list.k)>204800 ? ".gz" : "");
			writeFasta(output.resolve(name),list.words,list.k);
			int nr=0,ns=0,nb=0;
			for(int word:list.words){final boolean r=Arrays.binarySearch(list.real,word)>=0,s=Arrays.binarySearch(list.shadow,word)>=0;
				if(r){nr++;}if(s){ns++;}if(r && s){nb++;}}
			summary.append(name).append('\t').append(list.k).append('\t').append(list.words.length).append('\t').append(list.name)
				.append('\t').append(nr).append('\t').append(ns).append('\t').append(nb).append("\tunmeasured\n");
		}
		writeText("lists.tsv",summary.toString());writeText("recipe.txt",recipe.toString());
		verifyInputs();final String pass="PROTEIN_KMER_CANDIDATES_PASS lists="+lists.size()+" adopted=false";
		writeText("PASS",pass+"\n");outstream.println(pass);
	}

	/** Preserve accepted raw rankings without reading counts or selecting words again. */
	private void exportRanked() throws Exception{
		requireBgzip9();
		require(!text(rankedExport.resolve("accepted.utc")).trim().isEmpty(),"Missing candidate acceptance");
		for(String name:new String[]{"native.exit","exit"}){require(text(rankedExport.resolve(name)).equals("0\n"),"Failed candidate generation");}
		final Path result=rankedExport.resolve("result");
		final String[] manifest=text(result.resolve("lists.tsv")).split("\n",-1);
		require(manifest.length>=3 && manifest[0].equals(ProteinKmerWeightedCandidates.HEADER) && manifest[manifest.length-1].isEmpty(),"Wrong candidate manifest");
		final int total=manifest.length-2;
		require(text(result.resolve("PASS")).equals("PROTEIN_KMER_WEIGHTED_CANDIDATES_PASS lists="+total+" adopted=false\n"),"Wrong candidate result PASS");
		require(text(rankedExport.resolve("PASS")).equals("PROTEIN_KMER_CANDIDATES_RUN_PASS lists="+total+" adopted=false\n"),"Wrong candidate run PASS");
		final ArrayList<String> names=new ArrayList<String>(),payloads=new ArrayList<String>();
		final HashSet<String> seen=new HashSet<String>();
		final StringBuilder exported=new StringBuilder("file\tentries\tselection\torder\n");
		for(int row=1;row<=total;row++){
			final String[] f=manifest[row].split("\t",-1);require(f.length==11,"Malformed candidate manifest row");
			if(f[4].equals("reduced")){continue;}
			require(f[4].equals("raw") && f[1].equals("6") && f[5].equals("none") && (f[3].equals("shadow") || f[3].equals("weighted")),"Require explicit raw K6 candidate policy");
			final int size=Integer.parseInt(f[2]);require(size>0 && size<=16384,"Ranked export supports up to16384 retained words per list");
			final String stem="raw_"+f[3]+"_k6_"+size;
			require(seen.add(f[3]) && (f[0].equals(stem+".faa") || f[0].equals(stem+".faa.gz")),"Duplicate selection or unsafe candidate filename");
			final long nw=Long.parseLong(f[6]),pw=Long.parseLong(f[7]),nn=Long.parseLong(f[8]),pn=Long.parseLong(f[9]);
			final boolean weighted=f[3].equals("weighted");
			require(nn>0 && pn>0 && (weighted ? nw>0 && pw>0 : nw==1 && pw==0),"Invalid source weights/windows");
			final ProteinKmerWeightedCandidates.Weights weights=weighted ? new ProteinKmerWeightedCandidates.Weights(nn,pn,nw,pw) : null;
			final Path table=result.resolve(stem+".scores.tsv.gz");final StableInput stable=track(table,inputs);
			final ByteFile reader=ByteFile.makeByteFile(table.toString(),false);final LineParser1 parser=new LineParser1('\t');
			final int[] words=new int[size];int count=0,previousWord=-1;long previousScore=Long.MAX_VALUE;
			final StringBuilder fasta=new StringBuilder();
			try{
				require(Arrays.equals(reader.nextLine(),"rank\tpacked_word\tnegative_occurrences\tpositive_occurrences\tinteger_score".getBytes(StandardCharsets.US_ASCII)),"Wrong score-table schema");
				for(byte[] line=reader.nextLine();line!=null;line=reader.nextLine()){
					parser.set(line);require(parser.terms()==5 && count<size && parser.parseInt(0)==count+1,"Invalid rank/count");
					final int word=parser.parseInt(1);final long n=parser.parseLong(2),p=parser.parseLong(3),score=parser.parseLong(4);
					require(word>=0 && word<(1<<30) && n>=0 && n<=nn && p>=0 && p<=pn && score>0 &&
						score==(weighted ? weights.score(n,p) : n),"Invalid word/count/score");
					require(score<previousScore || (score==previousScore && word>previousWord),"Scores must descend with packed-word ties ascending");
					words[count++]=word;previousWord=word;previousScore=score;
					fasta.append(">rank_").append(count).append(" packed_word=").append(word).append('\n');
					for(int shift=25;shift>=0;shift-=5){final int code=(word>>>shift)&31;require(code<AA.length(),"Noncanonical packed residue");fasta.append(AA.charAt(code));}
					fasta.append('\n');
				}
			}finally{require(!reader.close(),"Score-table read failed");}
			require(count==size,"Incomplete score table");stable.verify();
			final int[] original=readFasta(result.resolve(f[0]),6,size,null,null,inputs);
			require(Arrays.equals(sortedUnique(words),sortedUnique(original)),"Ranked score/FASTA membership differs");
			names.add(stem+".ranked.faa.gz");payloads.add(fasta.toString());
			exported.append(names.get(names.size()-1)).append('\t').append(size).append('\t').append(f[3]).append("\tscore_descending_packed_word_ascending\n");
		}
		require(!names.isEmpty(),"No raw candidate lists to export");verifyInputs();Files.createDirectory(output);
		for(int i=0;i<names.size();i++){writeTextArtifact(output.resolve(names.get(i)),payloads.get(i));}
		writeText("lists.tsv",exported.toString());
		writeText("recipe.txt","source_package="+result+"\noperation=ranked_raw_export_only\ncounts_read=false\nreselection=false\ncompression=BGZF9\nfile_hashes_computed=0\nadopted=false\n");
		verifyInputs();final String pass="PROTEIN_KMER_RANKED_EXPORT_PASS lists="+names.size()+" adopted=false";
		writeText("PASS",pass+"\n");outstream.println(pass);
	}

	private int[] accepted(Path root,String population,int k,StringBuilder recipe) throws Exception{
		final Path p=root.resolve(population+"_k"+k),result=p.resolve("result");
		require(!text(p.resolve("accepted.utc")).trim().isEmpty(),"Missing owner acceptance");
		for(String name:new String[]{"native.exit","exit","watch.exit"}){require(text(p.resolve(name)).trim().equals("0"),"Failed count: "+p.resolve(name));}
		require(text(p.resolve("PASS")).equals("PROTEIN_KMER_COUNTS_PASS\n"),"Missing count PASS");
		require(text(p.resolve("native.time.tsv")).equals(text(p.resolve("accepted.time.tsv"))),"Changed accepted timing");
		final String sizes=text(p.resolve("accepted.output_sizes.tsv"));
		require(text(p.resolve("output_sizes.tsv")).equals(sizes),"Changed accepted sizes");
		bindAcceptedSizes(result,sizes,k==6 && maxK6>LIMIT);
		require(text(result.resolve("PASS")).equals("PROTEIN_KMER_FREQUENCY_PASS\n"),"Missing frequency PASS");
		final String inputRecipe=text(result.resolve("recipe.txt"));
		require(inputRecipe.contains("\npopulation="+population+"\n") && inputRecipe.contains("\nk="+k+"\n") &&
			inputRecipe.contains("\ntop=2048\n"),"Count population/K/top binding differs");
		recipe.append("input_").append(population).append("_k").append(k).append('=').append(result).append('\n');
		final Ranked top=readRankedData(result.resolve("top.tsv.gz"),result.resolve("top.faa.gz"),k,LIMIT,inputs);
		if(k!=6 || maxK6<=LIMIT){return top.words;}
		final String[] lines=text(result.resolve("summary.tsv")).split("\n",-1);
		require(lines.length==3 && lines[2].isEmpty() && lines[0].equals("population\tk\trecords\tresidues\tvalid_words\tdistinct_words\tselected_words"),"Invalid frequency summary");
		final String[] fields=lines[1].split("\t",-1);
		require(fields.length==7 && fields[0].equals(population) && fields[1].equals("6") && fields[6].equals("2048"),"Frequency summary binding differs");
		require(Long.parseLong(fields[2])>0 && Long.parseLong(fields[3])>0,"Invalid frequency population");
		final long windows=Long.parseLong(fields[4]),rows=Long.parseLong(fields[5]);
		final Ranked full=rankCounts(result.resolve("counts.tsv.gz"),6,maxK6,rows,windows,inputs);
		for(int i=0;i<LIMIT;i++){require(full.words[i]==top.words[i] && full.counts[i]==top.counts[i],"Full-count/top2048 parity failure at rank "+(i+1));}
		recipe.append("full_").append(population).append("_k6_rows=").append(rows).append("\nfull_")
			.append(population).append("_k6_windows=").append(windows).append("\nfull_").append(population).append("_top2048_parity=exact\n");
		return full.words;
	}

	/** Bind every consumed result file to the accepted stored-byte census. */
	void bindAcceptedSizes(Path result,String receipt,boolean full) throws Exception{
		final String[] names={"counts.tsv.gz","top.tsv.gz","top.faa.gz","summary.tsv","recipe.txt","PASS"};
		final long[] sizes=new long[names.length];Arrays.fill(sizes,-1);
		final String[] lines=receipt.split("\n",-1);
		require(lines.length==names.length+2 && lines[0].equals("file\tbytes") && lines[lines.length-1].isEmpty(),"Invalid accepted size receipt schema/row count");
		for(int row=1;row<lines.length-1;row++){
			final String[] fields=lines[row].split("\t",-1);
			require(fields.length==2 && fields[1].matches("[0-9]+"),"Invalid accepted size row");
			final int index=Arrays.asList(names).indexOf(fields[0]);
			require(index>=0 && sizes[index]<0,"Unknown or duplicate accepted size row: "+fields[0]);
			sizes[index]=Long.parseLong(fields[1]);require(sizes[index]>0,"Invalid accepted file size");
		}
		for(int i=0;i<names.length;i++){require(sizes[i]>=0,"Missing accepted size row: "+names[i]);}
		for(int i:(full ? new int[]{0,1,2,3,4,5} : new int[]{1,2,4,5})){
			final StableInput stable=track(result.resolve(names[i]),inputs);
			require(stable.before.size()==sizes[i],"Consumed file differs from accepted size: "+stable.path);
			stable.verify();
		}
	}

	/** Alternating encounter order; sorted result is emission order, not selection order. */
	public static int[] mixed(int[] real,int[] shadow,int limit){
		require(real!=null && shadow!=null && limit>0,"Invalid mixed selection");
		final int[] selected=new int[limit];int size=0;
		final IntHashSet seen=new IntHashSet(Math.min(limit,40960)*2);
		for(int rank=0;rank<Math.max(real.length,shadow.length);rank++){
			for(int side=0;side<2;side++){
				final int[] source=side==0 ? real : shadow;
				if(rank>=source.length){continue;}
				final int word=source[rank];
				if(!seen.add(word)){continue;}
				selected[size++]=word;
				if(size==limit){Arrays.sort(selected);return selected;}
			}
		}
		throw new IllegalArgumentException("Insufficient unique mixed candidates: "+size+" < "+limit);
	}

	public static int[] readRanked(Path table,Path fasta,int k,int expected) throws Exception{
		return readRankedData(table,fasta,k,expected,new ArrayList<StableInput>()).words;
	}
	static Ranked readRankedData(Path table,Path fasta,int k,int expected,ArrayList<StableInput> inputs) throws Exception{
		require((k==5 || k==6) && expected>0,"Invalid ranked list domain");
		final StableInput stable=track(table,inputs);
		final int[] words=new int[expected];final long[] counts=new long[expected];
		final ByteFile reader=ByteFile.makeByteFile(table.toString(),false);
		final LineParser1 parser=new LineParser1('\t');int rows=0;
		try{
			final byte[] header=reader.nextLine();
			require(header!=null && Arrays.equals(header,"rank\tpacked_word\tword\toccurrences".getBytes(StandardCharsets.US_ASCII)),"Wrong top-table schema");
			for(byte[] line=reader.nextLine();line!=null;line=reader.nextLine()){
				require(rows<expected,"Excess top rows");parser.set(line);
				require(parser.terms()==4 && parser.parseInt(0)==rows+1 && parser.length(2)==k,"Wrong rank/word length/schema");
				int packed=0;
				for(int i=0;i<k;i++){final int aa=AA.indexOf((char)parser.parseByte(2,i));require(aa>=0,"Noncanonical word");packed=(packed<<5)|aa;}
				final long count=parser.parseLong(3);
				require(parser.parseInt(1)==packed && count>0,"Invalid packed key/count");
				require(rows==0 || count<counts[rows-1] || (count==counts[rows-1] && packed>words[rows-1]),"Wrong frequency/tie order");
				words[rows]=packed;counts[rows++]=count;
			}
		}finally{require(!reader.close(),"Top-table read/close failed");}
		require(rows==expected,"Missing top rows");sortedUnique(words);stable.verify();
		readFasta(fasta,k,expected,words,counts,inputs);stable.verify();return new Ranked(words,counts);
	}


	/** Stream all rows while retaining at most limit primitive key/count pairs. */
	public static int[] rankCounts(Path table,int k,int limit,long expectedRows,long expectedWindows) throws Exception{
		return rankCounts(table,k,limit,expectedRows,expectedWindows,new ArrayList<StableInput>()).words;
	}
	private static Ranked rankCounts(Path table,int k,int limit,long expectedRows,long expectedWindows,ArrayList<StableInput> inputs) throws Exception{
		require((k==5 || k==6) && limit>0 && expectedRows>=limit && expectedRows<=ProteinKmerFrequency.universe(k) && expectedWindows>=expectedRows,"Invalid full-count selection domain");
		final StableInput stable=track(table,inputs);final TopHeap heap=new TopHeap(limit);
		final ByteFile reader=ByteFile.makeByteFile(table.toString(),false);final LineParser1 parser=new LineParser1('\t');
		long rows=0,windows=0;int previous=-1;
		try{
			final byte[] header=reader.nextLine();require(header!=null && Arrays.equals(header,"packed_word\tword\toccurrences".getBytes(StandardCharsets.US_ASCII)),"Wrong count-table schema");
			for(byte[] line=reader.nextLine();line!=null;line=reader.nextLine()){
				parser.set(line);require(parser.terms()==3 && parser.length(1)==k,"Malformed count row");
				int word=0;
				for(int i=0;i<k;i++){
					final int b=parser.parseByte(1,i),aa=b>=0 && b<ENCODE.length ? ENCODE[b] : -1;
					require(aa>=0,"Noncanonical count word");word=(word<<5)|aa;
				}
				final long count=parser.parseLong(2);
				require(word>previous && parser.parseInt(0)==word && count>0,"Count key/order/occurrences invalid");
				previous=word;rows++;windows=Math.addExact(windows,count);
				require(rows<=expectedRows && windows<=expectedWindows,"Count table exceeds accepted totals");heap.add(word,count);
			}
		}finally{require(!reader.close(),"Full count read/close failed");}
		require(rows==expectedRows && windows==expectedWindows,"Count table differs from accepted totals");
		stable.verify();return heap.ranked();
	}

	/** Worst-first primitive heap; larger keys lose equal-count ties. */
	static final class TopHeap{
		TopHeap(int capacity){words=new int[capacity];counts=new long[capacity];}
		void add(int word,long count){
			if(size<words.length){
				int at=size++;words[at]=word;counts[at]=count;
				while(at>0){final int parent=(at-1)/2;if(!worse(at,parent)){break;}swap(at,parent);at=parent;}
			}else if(count>counts[0] || (count==counts[0] && word<words[0])){words[0]=word;counts[0]=count;sink();}
		}
		Ranked ranked(){
			require(size==words.length,"Insufficient ranked words");final int[] keys=new int[size];final long[] values=new long[size];
			for(int i=size-1;i>=0;i--){keys[i]=words[0];values[i]=counts[0];--size;words[0]=words[size];counts[0]=counts[size];if(size>0){sink();}}
			return new Ranked(keys,values);
		}
		private void sink(){
			int parent=0;
			for(int child=1;child<size;child=parent*2+1){
				if(child+1<size && worse(child+1,child)){child++;}if(!worse(child,parent)){break;}swap(child,parent);parent=child;
			}
		}
		private boolean worse(int a,int b){return counts[a]<counts[b] || (counts[a]==counts[b] && words[a]>words[b]);}
		private void swap(int a,int b){final int word=words[a];words[a]=words[b];words[b]=word;final long count=counts[a];counts[a]=counts[b];counts[b]=count;}
		final int[] words;final long[] counts;private int size;
	}
	static final class Ranked{
		Ranked(int[] words_,long[] counts_){words=words_;counts=counts_;}
		final int[] words;final long[] counts;
	}
	private static int[] encoding(){
		final int[] codes=new int[128];Arrays.fill(codes,-1);for(int i=0;i<AA.length();i++){codes[AA.charAt(i)]=i;}return codes;
	}
	static long fastaBytes(int records,int k){
		long bytes=(long)records*(k+3);for(int i=1;i<=records;i++){for(int value=i;value>0;value/=10){bytes++;}}return bytes;
	}

	private static int[] readFasta(Path fasta,int k,int expected,int[] words,long[] counts,ArrayList<StableInput> inputs) throws Exception{
		final StableInput stable=track(fasta,inputs);Shared.AMINO_IN=true;
		final FileFormat format=FileFormat.testInput(fasta.toString(),FileFormat.FASTA,null,false,false);
		require(format.fasta() && format.amino(),"Require protein FASTA");
		final Streamer reader=StreamerFactory.makeStreamer(format,null,0,true,-1,false,true,1);
		final int[] result=new int[expected];int records=0;boolean terminal=false;reader.start();
		try{
			for(ListNum<Read> batch=reader.nextList();batch!=null;batch=reader.nextList()){
				for(Read read:batch){
					require(records<expected && read.mate==null && read.id!=null && !read.id.isEmpty() && read.bases!=null && read.bases.length==k,"Malformed FASTA population");
					int packed=0;
					for(byte b:read.bases){final int aa=AA.indexOf((char)b);require(aa>=0,"Noncanonical FASTA word");packed=(packed<<5)|aa;}
					if(words!=null){require(packed==words[records] && read.id.equals((records+1)+"_count_"+counts[records]),"FASTA/top-table mismatch");}
					result[records++]=packed;
				}
			}
			terminal=true;
		}finally{
			// FastaStreamer2ST.close is not cancellation: closing before its terminal races bf.nextLine.
			// Preserve the original validation exception while letting the producer finish normally.
			try{if(!terminal){while(reader.nextList()!=null){/* Drain bounded candidate FASTA input. */}}}
			finally{reader.close();}
		}
		require(!reader.errorState() && records==expected,"FASTA read/record count failure");
		sortedUnique(result);stable.verify();return result;
	}

	private static int[] sortedUnique(int[] words){
		final int[] sorted=words.clone();Arrays.sort(sorted);
		for(int i=1;i<sorted.length;i++){require(sorted[i]>sorted[i-1],"Duplicate input words");}return sorted;
	}
	private static int[] union(int[] a,int[] b){
		final int[] out=new int[a.length+b.length];int i=0,j=0,n=0;
		while(i<a.length || j<b.length){
			final int value=j==b.length || (i<a.length && a[i]<b[j]) ? a[i++] : b[j++];
			if(n==0 || value!=out[n-1]){out[n++]=value;}
		}return Arrays.copyOf(out,n);
	}
	static void writeFasta(Path path,int[] words,int k){
		writeFasta(path,words,k,null);
	}
	static void writeFasta(Path path, int[] words, int k, String header){
		writeArtifact(path, words, k, header, null);
	}
	static void writeTextArtifact(Path path,String text){
		writeArtifact(path, null, 0, null, text);
	}
	private static void writeArtifact(Path path, int[] words, int k, String header, String text){
		final boolean compressed=path.toString().endsWith(".gz");
		final boolean bgzf=ReadWrite.USE_BGZF,bgzip=ReadWrite.USE_BGZIP,force=ReadWrite.FORCE_BGZIP,
			nativeBgzf=ReadWrite.ALLOW_NATIVE_BGZF,change=ReadWrite.ALLOW_ZIPLEVEL_CHANGE;
		final int level=ReadWrite.ZIPLEVEL;
		try{
			if(compressed){
				requireBgzip9();
				ReadWrite.USE_BGZF=ReadWrite.USE_BGZIP=ReadWrite.FORCE_BGZIP=true;
				ReadWrite.ALLOW_NATIVE_BGZF=ReadWrite.ALLOW_ZIPLEVEL_CHANGE=false;ReadWrite.ZIPLEVEL=9;
			}
			final ByteStreamWriter writer=new ByteStreamWriter(FileFormat.testOutput(path.toString(),text==null ? FileFormat.FASTA : FileFormat.TEXT,null,compressed,false,false,false));
			final ByteBuilder row=new ByteBuilder();writer.start();
			try{
				if(text!=null){writer.print(text);}
				else{for(int i=0;i<words.length;i++){
					row.clear().append('>').append(i==0 && header!=null ? header : Integer.toString(i+1)).nl();
					for(int shift=(k-1)*5;shift>=0;shift-=5){row.append(AA.charAt((words[i]>>>shift)&31));}
					writer.print(row.nl());
				}}
			}finally{require(!writer.poisonAndWait(),"Candidate FASTA write failed");}
		}finally{
			ReadWrite.USE_BGZF=bgzf;ReadWrite.USE_BGZIP=bgzip;ReadWrite.FORCE_BGZIP=force;
			ReadWrite.ALLOW_NATIVE_BGZF=nativeBgzf;ReadWrite.ALLOW_ZIPLEVEL_CHANGE=change;ReadWrite.ZIPLEVEL=level;
		}
	}
	static void requireBgzip9(){
		require(Data.BGZIP() && Data.BGZIP_VERSION_levelFlag && Data.BGZIP_VERSION_threadsFlag,
			"Extended candidates require external bgzip supporting explicit -l 9 and -@ 1; native BGZF is not a level-9 fallback");
	}
	String text(Path path) throws Exception{
		final StableInput stable=track(path,inputs);
		require(stable.before.size()<=1048576,"Unexpectedly large receipt: "+path);
		final String text=new String(Files.readAllBytes(path),StandardCharsets.UTF_8);stable.verify();return text;
	}
	void writeText(String name,String text) throws Exception{
		Files.write(output.resolve(name),text.getBytes(StandardCharsets.UTF_8),StandardOpenOption.CREATE_NEW,StandardOpenOption.WRITE);
	}
	private static StableInput track(Path path,ArrayList<StableInput> inputs) throws Exception{
		final StableInput stable=new StableInput(path);inputs.add(stable);return stable;
	}
	void verifyInputs() throws Exception{for(StableInput input:inputs){input.verify();}}

	/** Structural identity only: no digest is instantiated. Detect final symlinks and target changes. */
	public static final class StableInput{
		public StableInput(Path path_) throws Exception{
			path=path_.toAbsolutePath().normalize();before=attributes(path);target=path.toRealPath();verify();
		}
		public void verify() throws Exception{
			final BasicFileAttributes after=attributes(path);
			require(target.equals(path.toRealPath()) && before.size()==after.size() && before.lastModifiedTime().equals(after.lastModifiedTime()) &&
				Objects.equals(before.fileKey(),after.fileKey()),"Input changed during selection: "+path);
		}
		private static BasicFileAttributes attributes(Path path) throws Exception{
			final BasicFileAttributes a=Files.readAttributes(path,BasicFileAttributes.class,LinkOption.NOFOLLOW_LINKS);
			require(a.isRegularFile(),"Require regular nonsymlink input: "+path);return a;
		}
		private final Path path,target;private final BasicFileAttributes before;
	}
	private static final class Selection{
		Selection(String name_,int k_,int[] words_,int[] real_,int[] shadow_){name=name_;k=k_;words=words_;real=real_;shadow=shadow_;}
		final String name;final int k;final int[] words,real,shadow;
	}
	private static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	/** Preserve the original2048 controls; all other requested sizes are ordinary parameters. */
	public static int[] parseSizes(String value){
		require(value!=null && !value.isEmpty(),"Require a comma-separated k6sizes list");
		final String[] fields=value.split(",",-1);final int[] sizes=new int[fields.length+1];sizes[0]=LIMIT;
		final HashSet<Integer> seen=new HashSet<Integer>();
		for(int i=0;i<fields.length;i++){
			require(fields[i].matches("[0-9]+"),"Invalid K6 size: "+fields[i]);
			final int size=Integer.parseInt(fields[i]);
			require(size>0 && size<=ProteinKmerFrequency.universe(6) && seen.add(size),"Invalid or duplicate K6 size: "+size);
			sizes[i+1]=size;
		}
		Arrays.sort(sizes);int n=0;
		for(int size:sizes){if(n==0 || size!=sizes[n-1]){sizes[n++]=size;}}
		return Arrays.copyOf(sizes,n);
	}
	private String sizesText(){
		final StringBuilder s=new StringBuilder();for(int size:k6sizes){if(s.length()>0){s.append(',');}s.append(size);}return s.toString();
	}
	private final PrintStream outstream;
	private Path realRoot,shadowRoot,legacy,rankedExport,output;
	private String arms,alphabet,blacklistProvenance="protein_kmer_candidates";
	private long negativeWeight=4,positiveWeight=1;
	private boolean versionedBlacklistHeaders=false;
	private int[] k6sizes={2048};
	private int maxK6=2048,candidateK=6;
	private final ArrayList<StableInput> inputs=new ArrayList<StableInput>();
	private static final int LIMIT=2048;
	private static final String AA="ARNDCEQGHILKMFPSTWYV";
	private static final int[] ENCODE=encoding();
}
