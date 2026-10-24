package prok;

import java.util.Locale;
import java.util.ArrayList;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import map.ObjectMap;
import parse.LineParser1;
import parse.PreParser;
import prok.NcrnaClippingComparison.Contig;
import prok.NcrnaClippingComparison.Feature;
import prok.NcrnaClippingComparison.Genome;
import shared.Shared;
import structures.ByteBuilder;

/** Grades saved, assembly-scoped calls against reviewed fresh prok CM loci.
 * Reuses the existing matcher, duplicate/mask ordering and endpoint accounting.
 * Source manifests are rehashed by the committed producer, not by this reader.
 * @author Brian, Raiden */
public final class NcrnaFreshPanelGrade {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args,NcrnaFreshPanelGrade.class,false);
		String inputs=null,sources=null,contigs=null,contigSources=null,verified=null,models=null,truth=null,history=null,masks=null,runs=null,batchRoot=null,batchGff=null,out=null,srpUnion=null,srpOld=null;
		int expected=-1,expectedModels=-1,expectedOldSrp=-1,expectedSrp=-1;long expectedLoci=-1,expectedHistory=-1;
		for(String arg:pp.args){final int eq=arg.indexOf('=');require(eq>0,"Use flag=value");final String k=arg.substring(0,eq).toLowerCase(Locale.ROOT),v=arg.substring(eq+1);
			if(k.equals("inputs")){inputs=v;}else if(k.equals("sources")){sources=v;}else if(k.equals("contigs")){contigs=v;}else if(k.equals("contigsources")){contigSources=v;}
			else if(k.equals("verified")){verified=v;}else if(k.equals("models")){models=v;}else if(k.equals("truth")){truth=v;}else if(k.equals("history")){history=v;}
			else if(k.equals("masks")){masks=v;}else if(k.equals("runs")){runs=v;}else if(k.equals("outprefix")){out=v;}
			else if(k.equals("batchroot")){batchRoot=v;}
			else if(k.equals("batchgff")){batchGff=v;}
			else if(k.equals("expected")){expected=Integer.parseInt(v);}else if(k.equals("expectedmodels")){expectedModels=Integer.parseInt(v);}
			else if(k.equals("expectedloci")){expectedLoci=Long.parseLong(v);}else if(k.equals("expectedhistoryrows")){expectedHistory=Long.parseLong(v);}
			else if(k.equals("srpunion")){srpUnion=v;}else if(k.equals("srpold")){srpOld=v;}
			else if(k.equals("expectedoldsrp")){expectedOldSrp=Integer.parseInt(v);}else if(k.equals("expectedsrp")){expectedSrp=Integer.parseInt(v);}
			else{throw new IllegalArgumentException("Unknown flag: "+k);}}
		require(inputs!=null && sources!=null && contigs!=null && contigSources!=null && verified!=null && models!=null && truth!=null && history!=null && masks!=null && (runs!=null ^ batchRoot!=null) && out!=null
			&& expected>0 && expectedModels>0 && expectedLoci>0 && expectedHistory>=0,"Explicit prok schemas, reviewed truth, saved runs and denominators required");
		final String resolvedBatchGff=batchCallsPath(batchRoot, batchGff);
		Shared.setThreads(1);ReadWrite.setZipThreads(1);
		final NcrnaFreshPanelGrade tool=new NcrnaFreshPanelGrade();tool.binding.prokInputs=true;tool.binding.verified=EukPanelGffValidation.hashes(verified);
		for(String path:new String[]{inputs,sources,contigs,contigSources,models,truth,history,masks,runs}){if(path!=null){tool.binding.bound(path);}}
		tool.batchRoot=batchRoot;tool.batchGff=resolvedBatchGff;
		tool.binding.panel.hashes=EukPanelGffValidation.hashes(sources);tool.binding.panel.contigHashes=EukPanelGffValidation.hashes(contigSources);
		tool.binding.loadProkInputs(inputs,expected);tool.binding.loadProkContigs(contigs);tool.loadGenomes();tool.loadModels(models,expectedModels);
		tool.grade.loadFrozenTruth(truth,expectedLoci,EukPanelTruthDisposition.PROK_POLICY);
		if(srpUnion!=null){
			require(srpOld!=null && expectedOldSrp>0 && expectedSrp>0,"SRP override requires its original IDs and old/new denominators");
			tool.binding.bound(srpUnion);tool.binding.bound(srpOld);NcrnaSrpUnionTruth.apply(tool.grade,srpUnion,srpOld,expectedOldSrp,expectedSrp);tool.unionMode=true;
		}else{require(srpOld==null && expectedOldSrp<0 && expectedSrp<0,"SRP override arguments must be supplied together");}
		tool.grade.loadHistoricalMasks(history,expectedHistory);tool.grade.loadMasks(masks);
		tool.run(runs,out,expected);Shared.closeStream(pp.outstream);
	}
	void loadGenomes(){
		assert(grade.genomes.isEmpty() && !binding.panel.assemblies.isEmpty()):"One validated prok universe supplies assembly-scoped contigs";
		for(EukPanelGffValidation.Assembly a:binding.panel.assemblies){
			final String fasta=binding.sourceFastas.get(a.id);final Genome g=new Genome(a.id,a.group,fasta,binding.bound(fasta));grade.addGenome(g);
			for(EukPanelGffValidation.Contig source:a.literal.values()){if(source!=null){
				require(source.length<=Integer.MAX_VALUE,"Existing interval matcher requires int-sized contigs");
				final Contig c=new Contig(g,source.id,(int)source.length,source.compartment);require(g.contigs.put(c.id,c)==null,"Duplicate native contig");g.ordered.add(c);
			}}
			// Native hash-map iteration must not determine match/audit order.
			g.ordered.sort((a0,b0)->a0.id.compareTo(b0.id));
		}
	}
	void loadModels(String path,int expected){
		assert(grade.modelFamilies.isEmpty()):"Model accessions define truth folds, never caller outcomes";
		final EukPanelLoci roster=new EukPanelLoci();roster.loadModels(path,expected);
		for(EukPanelLoci.Model model:roster.models.values()){if(model!=null){
			require(NcrnaClippingComparison.known(model.family) && grade.modelFamilies.put(model.accession,model.family)==null,"Duplicate/unknown model fold");
		}}
	}
	void run(String path,String prefix,int expected){
		assert(grade.frozenMode && grade.historyLoaded):"Matching requires reviewed fresh truth and retained historical masks";
		RrnaResourceIO.fresh(prefix,".summary.tsv",".strata.tsv",".audit.tsv.gz",".unknown.tsv.gz",".runtime.tsv",".workload.tsv",".call_overlaps.tsv.gz");
		if(unionMode){RrnaResourceIO.fresh(prefix,".srp_agreement.tsv");}
		ByteFile in=null;final ByteStreamWriter[] writers=new ByteStreamWriter[6];int count=0;long calls=0,srpMatched=0,srpAgree=0;
		final LineParser1 p=new LineParser1('\t');final ByteBuilder b=new ByteBuilder();
		try{
			final ByteStreamWriter audit=writers[0]=RrnaResourceIO.writer(prefix+".audit.tsv.gz"),unknown=writers[1]=RrnaResourceIO.writer(prefix+".unknown.tsv.gz");
			final ByteStreamWriter time=writers[2]=RrnaResourceIO.writer(prefix+".runtime.tsv"),work=writers[3]=RrnaResourceIO.writer(prefix+".workload.tsv"),overlaps=writers[4]=RrnaResourceIO.writer(prefix+".call_overlaps.tsv.gz");
			audit.println("arm\tpanel\tassembly\tcontig\tcompartment\tview\trole\tid\tlabel\tstart\tstop\tstrand\toutcome\tpartner\td5\td3\tgeometry\treference_kind");
			unknown.println("arm\tpanel\tassembly\tcontig\trole\tid\tstart\tstop\tstrand\treason\toriginal_attributes");
			time.println("arm\tpanel\ttask\tgenomes\tprocess_wall_seconds\tuser_seconds\tsystem_seconds");
			work.println("arm\tpanel\ttask\tmetric\tfamily\tcount");overlaps.println("arm\tpanel\tassembly\tcontig\tview\tcall_a\tcall_b\toverlap_bp\texact_coordinates\tstatus");
			if(unionMode){writers[5]=RrnaResourceIO.writer(prefix+".srp_agreement.tsv");writers[5].println("domain\tassembly\tcontig\ttruth_locus\ttruth_subtype\tcall_id\tcall_subtype\tagrees");}
			if(batchRoot==null){in=RrnaResourceIO.open(path);RrnaResourceIO.header(in,RUN_HEADER,path);}
			else{loadBatchCalls(unknown,b);}
			for(int index=0; index<expected; index++){
				final EukPanelGffValidation.Assembly a=binding.panel.assemblies.get(index);
				final Genome g=grade.genomes.get(a.id);String root=null;
				if(batchRoot==null){
					final byte[] line=in.nextLine();require(line!=null,"Every selected assembly needs a saved caller row");p.set(line);
					require(p.terms()==6 && p.parseInt(0)==index+1 && p.termEquals(a.id,1),"Saved caller roster must retain exact selected index and assembly order");
					root=p.parseString(2);final String file=p.parseString(4);
					require(file.equals(root+"/calls.gff"),"Each saved GFF must retain its original per-assembly run root");
					binding.bound(root+"/MANIFEST.sha80",p.parseString(3));binding.bound(file,p.parseString(5));binding.bound(root+"/INPUTS.sha80");
					binding.bound(root+"/time.txt");binding.bound(root+"/stderr");
					final ObjectMap<String,String> inputs=EukPanelGffValidation.hashes(root+"/INPUTS.sha80");
					require(g.hash.equals(inputs.get(g.fasta)),"Saved call's source FASTA must equal its own validated assembly");
					for(String input:inputs.keys()){if(input!=null){binding.bound(input,inputs.get(input));}}
					grade.loadCalls(file,"SHIPPED",g.contigs,unknown,b);
				}
				count++;
				for(Contig c:g.ordered){
					for(Feature f:c.calls){if(f.family.equals("6S")){require(NcrnaClippingComparison.inView(f,"6S_RF00013") ^ NcrnaClippingComparison.inView(f,"6S_RF01685"),"Each emitted6S call needs an explicit subtype label");}calls++;}
					for(String view:VIEWS){grade.grade(c,"SHIPPED",view,audit,b);NcrnaClippingComparison.callOverlaps(c,"SHIPPED",view,overlaps,b);}
					if(unionMode){final long[] agreement=srpAgreement(c,writers[5],b);srpMatched+=agreement[0];srpAgree+=agreement[1];}c.calls.clear();
				}
				if(batchRoot==null){
					NcrnaClippingComparison.readTime(root+"/time.txt",time,b,"SHIPPED",g.panel,Integer.toString(count),1);
					NcrnaClippingComparison.readWorkload(root+"/stderr",work,b,"SHIPPED",g.panel,Integer.toString(count),false);
				}
			}
			if(batchRoot==null){require(in.nextLine()==null,"Caller roster cannot contain extra assemblies beyond the frozen panel");}
			else{
				NcrnaClippingComparison.readTime(batchRoot+"/time.txt",time,b,"SHIPPED","ALL","batch",count);
				NcrnaClippingComparison.readWorkload(batchRoot+"/stderr",work,b,"SHIPPED","ALL","batch",false);
			}
			require(count==expected && count==grade.genomes.size(),"Every selected assembly must be graded exactly once");
			long accounted=0;for(String key:grade.totals.keys()){if(key!=null){final String view=key.substring(key.lastIndexOf('\t')+1);
				if(!view.equals("SRP") && !NcrnaClippingComparison.sixsView(view)){final long[] v=grade.totals.get(key);accounted+=v[0]+v[1]+v[7];}}}
			require(accounted==calls,"Disjoint family outcomes must conserve every emitted RNA call; union/subtype views are not additive");
		}finally{
			boolean error=in!=null && in.close();for(ByteStreamWriter w:writers){if(w!=null){error|=w.poisonAndWait();}}require(!error,"Fresh-panel input/output failed");
		}
		writeSummary(prefix);
		if(unionMode){System.err.println("SRP_SUBTYPE_AGREEMENT matched="+srpMatched+" agree="+srpAgree+" disagree="+(srpMatched-srpAgree));}
		System.err.println("FRESH_PROK_GRADE_COMPLETE assemblies="+count+" loci="+grade.referenceRows+" calls="+calls+" scope=INDEPENDENT-WITH-STATED-GAP no_caller_runs=true");
	}
	/** One saved combined call stream, with every source bound and contig identity
	 * unambiguous across the complete frozen panel. No per-genome timing is invented. */
	private void loadBatchCalls(final ByteStreamWriter unknown, final ByteBuilder b){
		require(batchRoot!=null && !grade.genomes.isEmpty(),"Batch grading needs an explicit saved run and loaded frozen assemblies");
		for(String file:new String[]{"MANIFEST.sha80","INPUTS.sha80","inputs.tsv","time.txt","stderr"}){binding.bound(batchRoot+'/'+file);}
		binding.bound(batchGff);
		final ByteFile roster=RrnaResourceIO.open(batchRoot+"/inputs.tsv");final LineParser1 p=new LineParser1('\t');int count=0;
		try{
			RrnaResourceIO.header(roster,"index\tassembly\tfasta\tsha80",batchRoot);
			for(byte[] row;(row=roster.nextLine())!=null;){
				p.set(row);require(p.terms()==4 && count<binding.panel.assemblies.size() && p.parseInt(0)==count+1,
					"Batch roster must retain every selected assembly exactly once in index order");
				final EukPanelGffValidation.Assembly a=binding.panel.assemblies.get(count++);final Genome g=grade.genomes.get(a.id);
				require(p.termEquals(g.id,1) && p.termEquals(g.fasta,2) && p.termEquals(g.hash,3),
					"Batch input identity must match its frozen assembly and FASTA: "+a.id);
			}
			require(count==binding.panel.assemblies.size(),"Batch roster must cover the complete frozen fresh panel");
		}finally{require(!roster.close(),"Batch input roster read failed");}
		final ObjectMap<String,String> inputs=EukPanelGffValidation.hashes(batchRoot+"/INPUTS.sha80");
		for(String input:inputs.keys()){if(input!=null){binding.bound(input,inputs.get(input));}}
		final ObjectMap<String,Contig> byContig=new ObjectMap<String,Contig>(String.class,Contig.class);
		for(Genome g:grade.ordered){
			require(g.hash.equals(inputs.get(g.fasta)),"Batch source FASTA must equal its validated assembly: "+g.id);
			for(Contig c:g.ordered){require(c.calls.isEmpty() && byContig.put(c.id,c)==null,
				"Combined fresh-panel calls require unique contig IDs across assemblies: "+c.id);}
		}
		grade.loadCalls(batchGff, "SHIPPED", byContig, unknown, b);
	}
	/** Preserve the original saved-run path and compression suffix. Never guess
	 * another run or silently prefer one of two possible GFF files. */
	public static String batchCallsPath(String root, String override){
		if(root==null){require(override==null, "batchgff requires batchroot; per-assembly runs retain their own paths");return null;}
		final String path=override==null ? root+"/calls.gff" : override;
		require(path.equals(root+"/calls.gff") || path.equals(root+"/calls.gff.gz"), "Batch GFF must be its original calls.gff or calls.gff.gz under batchroot");
		return path;
	}
	/** Report subtype labels on the same deterministic physical-union matching. */
	public static long[] srpAgreement(Contig c,ByteStreamWriter out,ByteBuilder b){
		assert(c!=null && out!=null):"Subtype agreement is a separate report on existing physical-locus matches";
		final ArrayList<Feature> truth=new ArrayList<Feature>(),calls=new ArrayList<Feature>();
		final ArrayList<NcrnaGffFamilyGrader.Locus> tl=new ArrayList<NcrnaGffFamilyGrader.Locus>(),cl=new ArrayList<NcrnaGffFamilyGrader.Locus>();
		for(Feature f:c.truth){if(NcrnaClippingComparison.inView(f,"SRP")){truth.add(f);tl.add(f.locus("SRP"));}}
		for(Feature f:c.calls){if(NcrnaClippingComparison.inView(f,"SRP")){calls.add(f);cl.add(f.locus("SRP"));}}
		final NcrnaGffFamilyGrader.MatchResult matches=NcrnaGffFamilyGrader.matchOneToOne(cl,tl,.5f,c.circular?c.length:0);final long[] result=new long[2];
		for(int i=0;i<truth.size();i++){final int j=matches.truthMatchedCallIdx[i];if(j<0){continue;}final Feature t=truth.get(i),q=calls.get(j);final boolean same=t.family.equals(q.family);result[0]++;if(same){result[1]++;}
			out.print(b.clear().append(c.genome.panel).tab().append(c.genome.id).tab().append(c.id).tab().append(t.id).tab().append(t.family).tab().append(q.id).tab().append(q.family).tab().append(same?1:0).nl());}
		return result;
	}
	void writeSummary(String prefix){
		assert(!grade.totals.isEmpty()):"Publish only after complete saved-call grading and conservation checks";
		final ObjectMap<String,long[]> totals=new ObjectMap<String,long[]>(String.class,long[].class);
		for(String key:grade.totals.keys()){if(key!=null){final int first=key.indexOf('\t'),last=key.lastIndexOf('\t');final String domain=key.substring(first+1,key.indexOf('\t',first+1)),view=key.substring(last+1);
			add(totals,domain+'\t'+display(view),grade.totals.get(key));add(totals,"all\t"+display(view),grade.totals.get(key));}}
		ByteStreamWriter summary=null,strata=null;final ByteBuilder b=new ByteBuilder();
		try{
			summary=RrnaResourceIO.writer(prefix+".summary.tsv");strata=RrnaResourceIO.writer(prefix+".strata.tsv");
			summary.println("# CM-reference-relative outcomes. SRP and6S union/subtype rows overlap; do not sum them. Unresolved historical6S masks apply only to the union, not to an inferred subtype. "+(batchRoot==null ? "Per-genome JVM startup is included in runtime.tsv." : "One batch process is reported once in runtime.tsv; no per-genome timing is inferred."));
			summary.println("scope\tdomain\tfamily\tTP\tFP\tFN\texact5\texact3\tbothExact\tduplicateFP\tunscoredCalls\tboundaryEligibleTP\tbothEndsEligibleTP\tprecision\trecall\texact5Fraction\texact3Fraction\tbothExactFraction");
			for(String key:EukPanelCmProposal.sorted(totals)){final long[] v=totals.get(key);b.clear().append("INDEPENDENT-WITH-STATED-GAP").tab().append(key);for(long n:v){b.tab().append(n);}
				NcrnaClippingComparison.fraction(b,v[0],v[0]+v[1]);NcrnaClippingComparison.fraction(b,v[0],v[0]+v[2]);NcrnaClippingComparison.fraction(b,v[3],v[8]);NcrnaClippingComparison.fraction(b,v[4],v[8]);NcrnaClippingComparison.fraction(b,v[5],v[9]);summary.print(b.nl());}
			strata.println("arm\tdomain\tcompartment\tview\tgeometry\tTP\tFN\texact5\texact3\tbothExact\tboundaryEligibleTP\tbothEndsEligibleTP");
			for(String key:EukPanelCmProposal.sorted(grade.truthStrata)){b.clear().append(key);for(long n:grade.truthStrata.get(key)){b.tab().append(n);}strata.print(b.nl());}
		}finally{boolean error=false;if(summary!=null){error|=summary.poisonAndWait();}if(strata!=null){error|=strata.poisonAndWait();}require(!error,"Fresh-panel summary write failed");}
	}
	static void add(ObjectMap<String,long[]> out,String key,long[] values){
		assert(values.length==10):"Frozen matcher supplies ten conserved detection/boundary counters";
		long[] total=out.get(key);if(total==null){total=new long[10];out.put(key,total);}for(int i=0;i<10;i++){total[i]+=values[i];}
	}
	static String display(String view){return view.equals("SSU")?"16S":view.equals("LSU")?"23S":view;}
	static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	public static final String RUN_HEADER="index\tassembly\tresult_root\tresult_sha80\tcalls_gff\tcalls_sha80";
	//TODO: Euk-enabled fresh-panel calls need an explicit5.8S view before use:
	//loadCalls classifies r58, but this prok-only roster omits it and conservation
	//above correctly fails if one appears. Do not drop those calls to pass grading.
	public static final String[] VIEWS={"RNaseP","SRP_small","SRP_large","SRP","tmRNA","6S_RF00013","6S_RF01685","6S","SSU","LSU","5S","tRNA"};
	final EukPanelCmEvidence binding=new EukPanelCmEvidence();
	final NcrnaClippingComparison grade=new NcrnaClippingComparison();
	boolean unionMode;
	private String batchRoot, batchGff;
}
