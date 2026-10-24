package prok;

import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Locale;
import java.util.Random;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import gff.GffLine;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import stream.Writer;
import stream.WriterFactory;
import structures.ByteBuilder;
import structures.IntList;
import structures.ListNum;

/** Parameterized genomic shreds with fixed sampled anchors across flank settings.
 * Reads indexed, uncompressed FASTA ranges only; no caller or alignment runs.
 * Random eukaryotic background excludes annotated rRNA plus a declared buffer.
 * Background calls still require adjudication before being called false positives.
 * This is a sampled evaluation instrument, never a seed/model construction corpus.
 * @author Raiden
 */
public final class Euk18sShredSet {

	public static void main(String[] args)throws IOException{
		final PreParser pp=new PreParser(args,Euk18sShredSet.class,false);
		final Euk18sShredSet tool=new Euk18sShredSet();String input=null,prefix=null;
		for(String arg:pp.args){
			final int split=arg.indexOf('=');require(split>0,"Expected flag=value");
			final String key=arg.substring(0,split).toLowerCase(Locale.ROOT),value=arg.substring(split+1);
			if(key.equals("in")){input=value;}else if(key.equals("outprefix")){prefix=value;}
			else if(key.equals("flanks")){final String[] values=value.split(",");tool.flanks=new int[values.length];for(int i=0;i<values.length;i++){tool.flanks[i]=Integer.parseInt(values[i]);}}
			else if(key.equals("cap")){tool.cap=Integer.parseInt(value);}else if(key.equals("background")){tool.background=Integer.parseInt(value);}
			else if(key.equals("prokbackground")){tool.prokBackground=Integer.parseInt(value);}
			else if(key.equals("backgroundlength")){tool.backgroundLength=Integer.parseInt(value);}
			else if(key.equals("maskbuffer")){tool.maskBuffer=Integer.parseInt(value);}
			else if(key.equals("frozenregions")){tool.frozenRegions=value;}else if(key.equals("referenceflank")){tool.referenceFlank=Integer.parseInt(value);require(tool.referenceFlank>=0,"Explicit referenceflank must be nonnegative");}
			else if(key.equals("requirecoverage")){tool.requireCoverage=parse.Parse.parseBoolean(value);}
			else if(key.equals("seed")){tool.seed=Long.parseLong(value);}else{throw new IllegalArgumentException("Unknown option: "+key);}
		}
		require(input!=null && prefix!=null,"Declare the pinned genome manifest and fresh output prefix");
		Shared.setThreads(1);tool.run(input,prefix);Shared.closeStream(pp.outstream);
	}

	public void run(String input,String prefix)throws IOException{
		require(nextId==0,"One new output generation per generator instance");
		require(cap>0 && background>=0 && prokBackground>=-1 && backgroundLength>0 && flanks.length>0 && flanks.length<=10,"Positive sampling cap/length and at most ten flank settings required");
		Arrays.sort(flanks);for(int i=0;i<flanks.length;i++){require(flanks[i]>=0 && (i==0 || flanks[i]>flanks[i-1]),"Flanks must be distinct nonnegative lengths");}
		final int maxFlank=flanks[flanks.length-1];require(maskBuffer>=0 && (frozenRegions!=null || maskBuffer>=maxFlank),"D1 sampling requires a buffer covering the largest flank; frozen anchors instead require explicit nonnegative exclusion checks");
		require((frozenRegions==null)==(referenceFlank<0),"Frozen regions and their reference flank must be supplied together");
		if(frozenRegions!=null){loadFrozen(frozenRegions);}
		RrnaResourceIO.fresh(prefix,".regions.tsv",".truth.tsv",".summary.tsv");
		for(int flank:flanks){RrnaResourceIO.fresh(prefix,".f"+flank+".fa");}
		final ArrayList<Genome> genomes=loadManifest(input);final Random random=new Random(seed);
		final ByteStreamWriter regions=RrnaResourceIO.writer(prefix+".regions.tsv"),truth=RrnaResourceIO.writer(prefix+".truth.tsv"),summary=RrnaResourceIO.writer(prefix+".summary.tsv");
		final ArrayList<Output> outputs=new ArrayList<Output>();
		try{
			regions.println(REGION_HEADER);truth.println(TRUTH_HEADER);
			summary.println("clade\tassembly\trole\tflank\tshreds\tbases\tgffRrnaRows\tgffRrnaAbsentContig\tseed\tmaskBuffer");
			for(int flank:flanks){outputs.add(new Output(prefix+".f"+flank+".fa",flank));}
			for(Genome g:genomes){
				loadFai(g);loadTruth(g);loadGff(g);for(Contig c:g.contigs){c.prepareMask(maskBuffer);}
				final ArrayList<Sample> samples=frozenRegions==null?select(g,random,maxFlank):selectFrozen(g);
				try(RandomAccessFile fasta=new RandomAccessFile(g.fasta,"r")){
					require(fasta.read()=='>',"Indexed extraction requires plain FASTA, not compressed data: "+g.fasta);
					for(Output out:outputs){
						final long[][] counts=new long[ROLES.length][2];
						for(Sample sample:samples){
							final Contig c=sample.contig;final int lo=left(sample.start,out.flank),hi=right(sample.stop,out.flank,c.length);
							if(sample.role.equals(EUK_BACKGROUND)){require(!c.overlaps(lo,hi),"Every background extraction must remain outside the buffered rRNA mask");}
							final String id="shred_"+sample.id;
							out.add(new Read(fetch(fasta,c,lo,hi),null,id,sample.id));
							writeRegion(regions,g,sample,out.flank,lo,hi,id);writeTruth(truth,g,sample,out.flank,lo,hi,id);
							final long[] n=counts[roleIndex(sample.role)];n[0]++;n[1]+=hi-lo+1;
						}
						final ByteBuilder b=new ByteBuilder();
						for(int i=0;i<ROLES.length;i++){final long[] n=counts[i];summary.print(b.clear().append(g.clade).tab().append(g.assembly).tab().append(ROLES[i]).tab().append(out.flank).tab().append(n[0]).tab().append(n[1]).tab().append(g.gffRows).tab().append(g.absentGffRows).tab().append(seed).tab().append(maskBuffer).nl());}
					}
				}
				System.err.println("SHRED_GENOME_PASS assembly="+g.assembly+" anchors="+samples.size()+" contigs="+g.contigs.size());
			}
			if(requireCoverage){
				for(String clade:Euk18sLongLocusCensus.CLADES){require(positiveClades.contains(clade),"D1 needs sampled positives in every eukaryotic clade; absent="+clade);}
				for(int i=4;i<ROLES.length;i++){require(roleCounts[i]>0,"D1 requires explicit background/decoy coverage; absent="+ROLES[i]);}
			}
			if(frozenRegions!=null){require(frozen.size()==genomes.size() && nextId==frozenCount,"Every frozen anchor and source assembly must be extracted exactly once");}
		}finally{
			boolean failed=false;for(Output out:outputs){failed|=out.close();}RrnaResourceIO.finish(regions,truth,summary);require(!failed,"Shred FASTA writer failed");
		}
		System.err.println("EUK18S_SHREDS_PASS genomes="+genomes.size()+" anchors="+nextId+" flankSettings="+flanks.length+" newAlignments=0");
		if(frozenRegions!=null){System.err.println("FROZEN_ANCHORS_PASS referenceFlank="+referenceFlank+" anchors="+frozenCount+" maskBuffer="+maskBuffer+" resampled=false");}
	}

	/** Read the original selection once. Changing flank length must not change
	 * random-background membership, feature identity, or globally assigned IDs. */
	void loadFrozen(String path){
		require(referenceFlank>=0 && frozen.isEmpty(),"One explicitly identified reference-flank selection is required");
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{RrnaResourceIO.header(in,REGION_HEADER,path);for(byte[] row;(row=in.nextLine())!=null;){p.set(row);require(p.terms()==13,"Frozen anchors must retain the original region schema");if(p.parseInt(1)!=referenceFlank){continue;}
			final String id=p.parseString(0);require(id.equals("shred_"+frozenCount),"Original sample IDs/order must remain contiguous and unchanged: "+id);roleIndex(p.parseString(4));
			final FrozenAnchor a=new FrozenAnchor(p.parseString(2),p.parseString(4),p.parseString(5),p.parseString(12),p.parseInt(6),p.parseInt(7),p.parseInt(8),p.parseInt(9),frozenCount++);
			require(a.lo>0 && a.hi>=a.lo && a.start>=a.lo && a.stop<=a.hi && a.stop>=a.start && p.parseInt(10)==a.start-a.lo && p.parseInt(11)==a.hi-a.stop,"Saved source and actual-flank coordinates must conserve the anchor span");
			final String assembly=p.parseString(3);FrozenSamples rows=frozen.get(assembly);if(rows==null){rows=new FrozenSamples();frozen.put(assembly,rows);}rows.list.add(a);
		}}finally{require(!in.close(),"Frozen region read failed");}require(frozenCount>0,"No frozen anchors for the declared reference flank");
	}
	ArrayList<Sample> selectFrozen(Genome g){
		final FrozenSamples rows=frozen.get(g.assembly);require(rows!=null,"Genome manifest must match the frozen source assemblies: "+g.assembly);
		final ObjectMap<String,Sample> candidates=new ObjectMap<String,Sample>(String.class,Sample.class);
		for(Sample s:g.candidates){require(candidates.put(s.role+'\t'+s.feature.id,s)==null,"Feature roles/IDs must identify unique source candidates");}
		final ArrayList<Sample> result=new ArrayList<Sample>();
		for(FrozenAnchor a:rows.list){require(a.clade.equals(g.clade) && a.id==nextId,"Frozen anchors must preserve original clade, manifest order and sample IDs");
			final Contig c=g.byName.get(a.seqid);require(c!=null,"Frozen source contig must exist under its exact original FASTA name");bounds(c,a.start,a.stop);
			require(a.lo==left(a.start,referenceFlank) && a.hi==right(a.stop,referenceFlank,c.length),"Saved reference bounds must match the indexed source contig and actual edge clipping");
			final Sample s;if(a.role.equals(EUK_BACKGROUND) || a.role.equals("prok_background")){
				require(a.role.equals("prok_background")==g.clade.equals("prok") && a.anchor.matches("random_[0-9]+"),"Frozen background role and original random anchor identity must remain explicit");
				s=new Sample(new Feature(a.anchor,c,a.start,a.stop,"+",false,"background"),a.role);
			}else{s=candidates.get(a.role+'\t'+a.anchor);require(s!=null && s.contig==c && s.start==a.start && s.stop==a.stop,"Frozen positive/decoy must match its unchanged original source feature");}
			s.id=nextId++;result.add(s);roleCounts[roleIndex(s.role)]++;if(s.role.startsWith("positive_")){positiveClades.add(g.clade);}
		}
		require(!result.isEmpty(),"A frozen source assembly cannot silently lose all anchors");return result;
	}

	static ArrayList<Genome> loadManifest(String path){
		final ArrayList<Genome> result=new ArrayList<Genome>();final ObjectSet<String> seen=new ObjectSet<String>(String.class);
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{RrnaResourceIO.header(in,MANIFEST_HEADER,path);for(byte[] row;(row=in.nextLine())!=null;){
			p.set(row);require(p.terms()==7,"Genome manifest needs clade/assembly/fasta/fai/gff/truth/expectedTruth");
			final Genome g=new Genome(p.parseString(0),p.parseString(1),p.parseString(2),p.parseString(3),p.parseString(4),p.parseString(5),p.parseLong(6));
			require(g.expectedTruth>=0 && (!g.clade.equals("prok") || g.expectedTruth==0),"Expected source truth cardinality must be explicit, nonnegative, and zero for prok controls");
			require((g.clade.equals("prok") || Euk18sLongLocusCensus.validClade(g.clade)) && seen.add(g.assembly),"One recognized clade and one exact source per assembly");
			require(g.clade.equals("prok")==g.truth.equals("."),"Only prokaryotic negatives omit the eukaryotic CM truth ledger");result.add(g);
		}}finally{require(!in.close(),"Genome manifest read failed");}require(!result.isEmpty(),"An empty manifest is not a shred set");return result;
	}
	public static void loadFai(Genome g){
		final ByteFile in=RrnaResourceIO.open(g.fai);final LineParser1 p=new LineParser1('\t');
		try{for(byte[] row;(row=in.nextLine())!=null;){
			p.set(row);require(p.terms()==5,"Plain FASTA FAI requires five columns");final long n=p.parseLong(1);
			require(n>0 && n<Integer.MAX_VALUE,"Per-contig coordinates must fit positive Java int bounds");
			final Contig c=new Contig(p.parseString(0),(int)n,p.parseLong(2),p.parseInt(3),p.parseInt(4));
			require(c.offset>0 && c.lineBases>0 && c.lineBytes>=c.lineBases && c.lineBytes-c.lineBases<=2,"FAI must describe fixed-width lines with LF/CRLF separators");
			require(g.byName.put(c.name,c)==null,"Duplicate FASTA contig name would make coordinate extraction ambiguous");g.contigs.add(c);
			final String bare=accession(c.name);require(g.byAccession.put(bare,c)==null,"Duplicate accession alias in FASTA: "+bare);
		}}finally{require(!in.close(),"FAI read failed");}require(!g.contigs.isEmpty(),"No FASTA contigs in declared index");
	}
	public static void loadTruth(Genome g){
		if(g.truth.equals(".")){return;}
		final ByteFile in=RrnaResourceIO.open(g.truth);final LineParser1 p=new LineParser1('\t');final ObjectSet<String> seen=new ObjectSet<String>(String.class);long observed=0;
		try{RrnaResourceIO.header(in,RrnaSpeciesKeys.TRUTH_HEADER,g.truth);for(byte[] row;(row=in.nextLine())!=null;){
			p.set(row);require(p.terms()==22,"Expected complete 22-column CM truth ledger");if(!p.termEquals(g.assembly,2)){continue;}
			require(p.termEquals(g.clade,1) && seen.add(p.parseString(0)),"Truth clade/identity must bind this assembly without duplicates");
			observed++;
			final Contig c=g.resolve(p.parseString(4));require(c!=null && c.length==p.parseInt(9),"Truth contig and length must match exact FASTA index");
			final int start=p.parseInt(5),stop=p.parseInt(6);bounds(c,start,stop);require(stop-start+1==p.parseInt(8),"CM genomic span must conserve inclusive length");
			final String strand=p.parseString(7),comp=p.parseString(20);require(strand.equals("+") || strand.equals("-"),"CM strand must be explicit");Euk18sM1Grade.validComp(comp);
			final boolean complete=flag(p,17),construction=flag(p,21);
			require(complete==(p.parseInt(12)==1 && p.parseInt(13)==1831 && p.termEquals("no",14) && stop-start+1>=916),"Admission must match the recorded full-model/untruncated/half-length rule");
			require(construction==(complete && comp.equals("no_keyword")),"Construction flag must bind admitted no-keyword population");
			final Feature f=new Feature(p.parseString(0),c,start,stop,strand,complete,comp);c.truth.add(f);c.mask.add(new int[]{start,stop});
			if(construction){g.candidates.add(new Sample(f,positiveRole(stop-start+1)));}
			else if(complete && (comp.equals("mitochondrial") || comp.equals("plastid"))){g.candidates.add(new Sample(f,comp.equals("mitochondrial")?"decoy_mito_ssu":"decoy_plastid_ssu"));}
		}}finally{require(!in.close(),"CM truth read failed");}require(observed==g.expectedTruth,"CM truth count differs from manifest expectation for "+g.assembly+": "+observed+" versus "+g.expectedTruth);
	}
	public static void loadGff(Genome g){
		final ByteFile in=RrnaResourceIO.open(g.gff);final ObjectSet<String> decoys=new ObjectSet<String>(String.class);final LineParser1 p=new LineParser1('\t');
		boolean assemblyBound=false;final boolean needsAssembly=g.assembly.startsWith("GCF_") || g.assembly.startsWith("GCA_");
		try{for(byte[] row;(row=in.nextLine())!=null;){
			if(row.length==0){continue;}
			if(row[0]=='#'){
				final String header=new String(row,java.nio.charset.StandardCharsets.US_ASCII);
				if(needsAssembly && header.startsWith("#!genome-build-accession ")){require(header.trim().equals("#!genome-build-accession NCBI_Assembly:"+g.assembly),"GFF assembly header must match the declared versioned assembly: "+g.assembly);assemblyBound=true;}
				continue;
			}
			p.set(row);require(p.terms()==9,"Source annotation must be nine-column GFF");
			if(!p.termEquals("rRNA",2) && !p.termEquals("rRNA_gene",2) && !p.termEquals("pseudogenic_rRNA",2)
				&& !rrnaBiotype(p)){continue;}
			final GffLine line=new GffLine(row);
			final String attr=line.attributes==null?"":line.attributes.toLowerCase(Locale.ROOT);
			g.gffRows++;final Contig c=g.resolve(line.seqid);if(c==null){g.absentGffRows++;continue;}
			bounds(c,line.start,line.stop);c.mask.add(new int[]{line.start,line.stop});
			String role=null;if(attr.contains("product=16s")){role="decoy_16s";}
			else if(attr.contains("product=28s") || attr.contains("product=25s") || attr.contains("product=26s")){role="decoy_28s";}
			if(role!=null){require(line.strand==0 || line.strand==1,"Decoy needs a declared strand");
				final String key=c.name+':'+line.start+':'+line.stop+':'+line.strand+':'+role;
				if(decoys.add(key)){g.candidates.add(new Sample(new Feature(key,c,line.start,line.stop,line.strand==0?"+":"-",true,"annotation"),role));}
			}
		}}finally{require(!in.close(),"rRNA annotation read failed");}
		require(!needsAssembly || assemblyBound,"Versioned assembly requires its original GFF assembly header, not a filename guess");
	}
	public ArrayList<Sample> select(Genome g,Random random,int maxFlank){
		final ArrayList<Sample> result=new ArrayList<Sample>();final ObjectMap<String,Samples> buckets=new ObjectMap<String,Samples>(String.class,Samples.class);
		final ArrayList<String> roles=new ArrayList<String>();
		for(Sample sample:g.candidates){Samples bucket=buckets.get(sample.role);if(bucket==null){bucket=new Samples();buckets.put(sample.role,bucket);roles.add(sample.role);}bucket.list.add(sample);}
		Collections.sort(roles);
		for(String role:roles){final ArrayList<Sample> list=buckets.get(role).list;Collections.shuffle(list,random);for(int i=0;i<Math.min(cap,list.size());i++){result.add(list.get(i));}}
		final ArrayList<StartSpan> spans=startSpans(g,backgroundLength,maxFlank);
		final long starts=spans.isEmpty()?0:spans.get(spans.size()-1).cumulative;
		final int requested=g.clade.equals("prok") && prokBackground>=0?prokBackground:background;
		require(requested<=starts,"Not enough distinct buffered-rRNA-excluded start positions for the requested background count");
		final ObjectSet<String> sampled=new ObjectSet<String>(String.class);int made=0;
		for(long attempts=0;made<requested && attempts<10000L*Math.max(1,requested);attempts++){
			final long position=bounded(random,starts);int low=0,high=spans.size();
			while(low<high){final int mid=(low+high)>>>1;if(spans.get(mid).cumulative<=position){low=mid+1;}else{high=mid;}}
			final StartSpan span=spans.get(low);final Contig chosen=span.contig;
			final int start=span.start+(int)(position-(low==0?0:spans.get(low-1).cumulative)),stop=start+backgroundLength-1;
			require(g.clade.equals("prok") || !chosen.overlaps(left(start,maxFlank),right(stop,maxFlank,chosen.length)),"Complement sampler must exclude every buffered rRNA interval");
			if(!sampled.add(chosen.name+':'+start)){continue;}
			result.add(new Sample(new Feature("random_"+made,chosen,start,stop,"+",false,"background"),g.clade.equals("prok")?"prok_background":EUK_BACKGROUND));made++;
		}
		require(made==requested,"Not enough rRNA-excluded start positions; report the gap rather than shrinking the requested denominator");
		for(Sample sample:result){sample.id=nextId++;roleCounts[roleIndex(sample.role)]++;if(sample.role.startsWith("positive_")){positiveClades.add(g.clade);}}
		require(!result.isEmpty(),"Source produced no requested shreds: "+g.assembly);return result;
	}
	public static ArrayList<StartSpan> startSpans(Genome g,int coreLength,int maxFlank){
		require(coreLength>0 && maxFlank>=0,"Start-space geometry requires positive cores and nonnegative flanks");
		final ArrayList<StartSpan> spans=new ArrayList<StartSpan>();long total=0;
		for(Contig c:g.contigs){
			final long last=(long)c.length-coreLength+1;if(last<1){continue;}long next=1;
			if(!g.clade.equals("prok")){for(int i=0;i<c.starts.size;i++){
				// A core starting here would overlap a buffered annotation after max-flank extraction.
				final long from=Math.max(1L,(long)c.starts.get(i)-coreLength+1-maxFlank),to=Math.min(last,(long)c.ends.get(i)+maxFlank);
				if(from>next){final long end=Math.min(last,from-1);if(end>=next){total=Math.addExact(total,end-next+1);spans.add(new StartSpan(c,(int)next,total));}}
				next=Math.max(next,to+1);if(next>last){break;}
			}}
			if(next<=last){total=Math.addExact(total,last-next+1);spans.add(new StartSpan(c,(int)next,total));}
		}return spans;
	}
	public static byte[] fetch(RandomAccessFile fasta,Contig c,int start,int stop)throws IOException{
		bounds(c,start,stop);final byte[] result=new byte[stop-start+1];int target=0,position=start-1;
		while(target<result.length){final int col=position%c.lineBases,n=Math.min(result.length-target,c.lineBases-col);
			final long offset=Math.addExact(c.offset,Math.addExact(Math.multiplyExact((long)(position/c.lineBases),c.lineBytes),col));
			fasta.seek(offset);fasta.readFully(result,target,n);position+=n;target+=n;
		}
		for(int i=0;i<result.length;i++){byte b=result[i];if(b>='a' && b<='z'){b-=32;}require("ACGTURYSWKMBDHVN".indexOf(b)>=0,"FAI range contains a header, line ending or invalid base; wrong/stale index: "+c.name);result[i]=b;}return result;
	}
	static void writeRegion(ByteStreamWriter out,Genome g,Sample s,int flank,int lo,int hi,String id){
		out.print(new ByteBuilder().append(id).tab().append(flank).tab().append(g.clade).tab().append(g.assembly).tab().append(s.role).tab().append(s.contig.name).tab().append(lo).tab().append(hi).tab().append(s.start).tab().append(s.stop).tab().append(s.start-lo).tab().append(hi-s.stop).tab().append(s.feature.id).nl());
	}
	static void writeTruth(ByteStreamWriter out,Genome g,Sample sample,int flank,int lo,int hi,String id){
		final ByteBuilder b=new ByteBuilder();for(Feature f:sample.contig.truth){
			if(f.stop<lo || f.start>hi){continue;}final boolean clipped=f.start<lo || f.stop>hi;
			out.print(b.clear().append(id).tab().append(flank).tab().append(f.id).tab().append(Math.max(f.start,lo)-lo+1).tab().append(Math.min(f.stop,hi)-lo+1).tab().append(f.strand).tab().append(f.complete?"complete":"partial").tab().append(f.compartment).tab().append(clipped?1:0).tab().append(f.start).tab().append(f.stop).nl());
		}
	}
	static long bounded(Random random,long bound){require(bound>0,"Random start universe cannot be empty");long bits,value;do{bits=random.nextLong()>>>1;value=bits%bound;}while(bits-value+(bound-1)<0);return value;}
	static boolean rrnaBiotype(LineParser1 p){
		p.setBounds(8);final byte[] bytes=p.line();final int end=p.b();
		for(int start=p.a();start<end;){
			int stop=start;while(stop<end && bytes[stop]!=';'){stop++;}
			for(String key:BIOTYPE_KEYS){final int value=start+key.length();
				if(value<=stop && matches(bytes,start,value,key)){
					for(String kind:RRNA_BIOTYPES){if(stop-value==kind.length() && matches(bytes,value,stop,kind)){return true;}}
				}
			}start=stop+1;
		}return false;
	}
	static boolean matches(byte[] bytes,int start,int stop,String lower){assert(stop-start==lower.length()):"Literal attribute match requires exact token bounds";for(int i=start;i<stop;i++){int c=bytes[i];if(c>='A' && c<='Z'){c+=32;}if(c!=lower.charAt(i-start)){return false;}}return true;}
	static int roleIndex(String role){for(int i=0;i<ROLES.length;i++){if(ROLES[i].equals(role)){return i;}}throw new IllegalArgumentException("Unknown shred role: "+role);}
	static String positiveRole(int n){return n<1600?"positive_short":n<=2200?"positive_ordinary":n<=2500?"positive_long":"positive_very_long";}
	static String accession(String name){if(name.startsWith("tid|")){final int p=name.indexOf('|',4);require(p>4 && p<name.length()-1,"Malformed native accession header");return name.substring(p+1);}return name;}
	static boolean flag(LineParser1 p,int n){require(p.termEquals("0",n) || p.termEquals("1",n),"Truth population flags must be Boolean");return p.termEquals("1",n);}
	static int left(int start,int pad){return (int)Math.max(1L,(long)start-pad);}
	static int right(int stop,int pad,int length){return (int)Math.min(length,(long)stop+pad);}
	static void bounds(Contig c,int start,int stop){require(start>0 && stop>=start && stop<=c.length,"Inclusive feature/extraction bounds must lie in the indexed contig: "+c.name);}
	static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}

	public static final class Contig {
		public Contig(String name,int length,long offset,int bases,int bytes){this.name=name;this.length=length;this.offset=offset;lineBases=bases;lineBytes=bytes;}
		public void prepareMask(int buffer){
			assert(starts.size==0):"Merge each source mask once before sampling";
			mask.sort(Comparator.comparingInt(a->a[0]));
			for(int[] interval:mask){final int lo=left(interval[0],buffer),hi=right(interval[1],buffer,length);
				if(starts.size>0 && (long)lo<=ends.get(ends.size-1)+1L){ends.set(ends.size-1,Math.max(ends.get(ends.size-1),hi));}
				else{starts.add(lo);ends.add(hi);}}
		}
		public boolean overlaps(int lo,int hi){int a=0,b=starts.size;while(a<b){final int mid=(a+b)>>>1;if(starts.get(mid)<=hi){a=mid+1;}else{b=mid;}}return a>0 && ends.get(a-1)>=lo;}
		final String name;final int length,lineBases,lineBytes;final long offset;
		public final ArrayList<int[]> mask=new ArrayList<int[]>();final IntList starts=new IntList(),ends=new IntList();final ArrayList<Feature> truth=new ArrayList<Feature>();
	}
	static final class Feature {
		Feature(String id,Contig c,int start,int stop,String strand,boolean complete,String compartment){this.id=id;contig=c;this.start=start;this.stop=stop;this.strand=strand;this.complete=complete;this.compartment=compartment;}
		final String id,strand,compartment;final Contig contig;final int start,stop;final boolean complete;
	}
	static final class Sample {
		Sample(Feature f,String role){feature=f;contig=f.contig;start=f.start;stop=f.stop;this.role=role;}
		final Feature feature;final Contig contig;final int start,stop;final String role;long id;
	}
	static final class Samples {final ArrayList<Sample> list=new ArrayList<Sample>();}
	static final class FrozenSamples {final ArrayList<FrozenAnchor> list=new ArrayList<FrozenAnchor>();}
	static final class FrozenAnchor {
		FrozenAnchor(String c,String r,String s,String a,int lo_,int hi_,int start_,int stop_,long id_){clade=c;role=r;seqid=s;anchor=a;lo=lo_;hi=hi_;start=start_;stop=stop_;id=id_;}
		final String clade,role,seqid,anchor;final int lo,hi,start,stop;final long id;
	}
	public static final class StartSpan {StartSpan(Contig c,int s,long n){contig=c;start=s;cumulative=n;}final Contig contig;public final int start;public final long cumulative;}
	public static final class Genome {
		public Genome(String c,String a,String f,String i,String g,String t,long expected){clade=c;assembly=a;fasta=f;fai=i;gff=g;truth=t;expectedTruth=expected;}
		Contig resolve(String id){final Contig exact=byName.get(id);return exact!=null?exact:byAccession.get(id);}
		public final String clade,assembly,fasta,fai,gff,truth;final long expectedTruth;long gffRows,absentGffRows;
		public final ArrayList<Contig> contigs=new ArrayList<Contig>();final ArrayList<Sample> candidates=new ArrayList<Sample>();
		final ObjectMap<String,Contig> byName=new ObjectMap<String,Contig>(String.class,Contig.class),byAccession=new ObjectMap<String,Contig>(String.class,Contig.class);
	}
	public static final class Output {
		public Output(String path,int flank){this.flank=flank;writer=WriterFactory.makeWriter(FileFormat.testOutput(path,FileFormat.FASTA,null,true,false,false,false),null,1,null,false);writer.start();}
		void add(Read r){batch.add(r);if(batch.size()==64){flush();}}
		void flush(){if(!batch.isEmpty()){writer.addReads(new ListNum<Read>(batch,batchId++));batch=new ArrayList<Read>(64);}}
		public boolean close(){flush();return writer.poisonAndWait();}
		final int flank;final Writer writer;ArrayList<Read> batch=new ArrayList<Read>(64);long batchId;
	}
	public static final String MANIFEST_HEADER="clade\tassembly\tfasta\tfai\tgff\ttruth\texpectedTruth";
	public static final String REGION_HEADER="shred\tflank\tclade\tassembly\trole\tsourceSeqid\tsourceStart1\tsourceStop1\tanchorStart1\tanchorStop1\tactualLeftFlank\tactualRightFlank\tanchorId";
	public static final String TRUTH_HEADER="shred\tflank\tlocus\tstart1\tstop1\tstrand\tcompleteness\tcompartment\tclipped\tsourceStart1\tsourceStop1";
	static final String EUK_BACKGROUND="rrna_excluded_euk_background";
	static final String[] ROLES={"positive_short","positive_ordinary","positive_long","positive_very_long","decoy_16s","decoy_28s","decoy_mito_ssu","decoy_plastid_ssu",EUK_BACKGROUND,"prok_background"};
	static final String[] BIOTYPE_KEYS={"gene_biotype=","transcript_biotype=","biotype=","ncrna_class="},RRNA_BIOTYPES={"rrna","rrna_pseudogene","pseudogenic_rrna"};
	final long[] roleCounts=new long[ROLES.length];final ObjectSet<String> positiveClades=new ObjectSet<String>(String.class);
	public int[] flanks={0,500,1000};public int cap=32,background=640,prokBackground=-1,backgroundLength=1843,maskBuffer=1000;public long seed=20261005L,nextId;public boolean requireCoverage=true;
	public String frozenRegions;public int referenceFlank=-1;long frozenCount;
	final ObjectMap<String,FrozenSamples> frozen=new ObjectMap<String,FrozenSamples>(String.class,FrozenSamples.class);
}
