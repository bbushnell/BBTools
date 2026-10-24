package prok;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import prok.NcrnaGffFamilyGrader.Locus;
import prok.NcrnaGffFamilyGrader.MatchedPair;
import shared.Shared;
import structures.ByteBuilder;

/** One-to-one 18S M1 grading and observed-stage attribution. All included CM
 * loci are matched together before completeness/compartment stratification.
 * No-keyword is an audit label, not proof of nuclear origin.
 * @author Raiden
 */
public final class Euk18sM1Grade {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args,Euk18sM1Grade.class,false);args=pp.args;
		String calls=null,truth=null,stages=null,prefix=null,arm=null,clade=null,assembly=null,contigCatalog=null;int models=-1,callerMinLen=-1;boolean negative=false,diagnostics=true;
		for(String arg:args){final int e=arg.indexOf('=');require(e>0,"Expected flag=value");final String k=arg.substring(0,e).toLowerCase(Locale.ROOT),v=arg.substring(e+1);
			if(k.equals("calls")){calls=v;}else if(k.equals("truth")){truth=v;}else if(k.equals("stages")){stages=v;}else if(k.equals("outprefix")){prefix=v;}
			else if(k.equals("arm")){arm=v;}else if(k.equals("clade")){clade=v;}else if(k.equals("assembly")){assembly=v;}
			else if(k.equals("contigs")){contigCatalog=v;}else if(k.equals("callerminlen")){callerMinLen=Parse.parseIntKMG(v);}
			else if(k.equals("diagnostics")){diagnostics=Parse.parseBoolean(v);}
			else if(k.equals("models")){models=Parse.parseIntKMG(v);}else if(k.equals("negative")){negative=Parse.parseBoolean(v);}else{throw new IllegalArgumentException(k);}}
		require(calls!=null && truth!=null && prefix!=null && arm!=null && clade!=null && assembly!=null,"Required calls/truth/outprefix/arm/clade/assembly");
		require(diagnostics ? stages!=null && models>0 : stages==null && models==-1 && contigCatalog==null && callerMinLen==-1 && !negative,"Diagnostic mode requires stages/models; diagnostics=f consumes only whole-genome calls and truth, without invented stage or defline evidence");
		require(negative==clade.equals("prok"),"Only the prok cohort is an explicit negative euk18S control");
		require((contigCatalog==null)==(callerMinLen<0) && callerMinLen!=0,"Optional sealed contig catalog must bind the actual positive caller minLen");
		final Euk18sM1Grade tool=new Euk18sM1Grade(assembly,negative);
		if(diagnostics){tool.run(calls,truth,stages,prefix,arm,clade,models,contigCatalog,callerMinLen);}
		else{tool.runWithoutDiagnostics(calls,truth,prefix,arm,clade);}
		Shared.closeStream(pp.outstream);
	}
	public Euk18sM1Grade(String a,boolean negative){require(a!=null&&!a.isEmpty(),"Every grade names one assembly");assembly=a;grade=new RrnaSeedBenchmarkGrade("euk18S",negative);}
	/** Saved whole-genome output has no seed stream. Match all loci together,
	 * retain truth labels, and keep calls on unlabelled contigs explicitly unknown. */
	public void runWithoutDiagnostics(String calls,String truth,String prefix,String arm,String clade){
		require(!grade.negative,"Whole-genome truth includes zero-locus assemblies without declaring biological negatives");
		RrnaCmWindowParity.fresh(prefix,new String[]{".summary.tsv",".loci.tsv.gz",".strata.tsv",".evidence.tsv"});
		loadTruth(truth);loadCalls(calls,false);grade.grade();grade.write(prefix,arm,clade,assembly);
		writeStrata(prefix,arm,clade,WHOLE_COMPS);
		final ByteStreamWriter out=RrnaSourceTargetPrep.writer(prefix+".evidence.tsv");
		try{out.println("assembly\tmodel_evidence\tstage_evidence\tcompartment_evidence\tunmatched_interpretation");
			out.print(new ByteBuilder().append(assembly).append("\tRF01960_ONLY\tNA\tTRUTH_KEYWORDS_ONLY_OTHER_CONTIGS_UNOBSERVED\tREFERENCE_RELATIVE_NOT_PROVEN_BIOLOGICAL_FP\n"));
		}finally{RrnaCmWindowParity.finish(out);}
		System.err.println("EUK18S_WHOLE_GRADE_PASS assembly="+assembly+" TP="+grade.tp+" unmatched="+grade.fp+" FN="+grade.fn+" model_evidence=RF01960_ONLY stage_evidence=NA");
	}
	public void run(String calls,String truth,String stages,String prefix,String arm,String clade,int models){
		run(calls,truth,stages,prefix,arm,clade,models,null,-1);
	}
	public void run(String calls,String truth,String stages,String prefix,String arm,String clade,int models,String catalog,int callerMinLen){
		loadTruth(truth);loadCompartments(stages,catalog,callerMinLen);writeCompartmentAudit(prefix);loadCalls(calls);grade.grade();grade.write(prefix,arm,clade,assembly);
		writeStrata(prefix,arm,clade);
		final Euk5sStageFailureSummary failure=new Euk5sStageFailureSummary(assembly,models,"euk18S");
		require(grade.fn<=Integer.MAX_VALUE,"Stage reducer uses integer per-assembly FN cardinality");
		failure.loadLoci(prefix+".loci.tsv.gz",(int)grade.fn);failure.loadStages(stages);failure.write(prefix+".stage");
		writeFailures(prefix);
		System.err.println("EUK18S_M1_GRADE_PASS assembly="+assembly+" TP="+grade.tp+" FP="+grade.fp+" FN="+grade.fn);
	}
	public void loadTruth(String path){
		final ByteFile in=RrnaSourceTargetPrep.open(path);final LineParser1 p=new LineParser1('\t');
		try{header(in,TRUTH_HEADER,path);for(byte[] row;(row=in.nextLine())!=null;){p.set(row);
			require(p.terms()==8&&p.termEquals(assembly,0),"Truth must retain exact assembly identity and eight columns");
			final Locus l=new Locus(p.parseString(1),"euk18S",p.parseInt(2),p.parseInt(3),RrnaSeedBenchmarkGrade.strand(p.parseString(4)));
			final String kind=p.parseString(5),comp=p.parseString(6);require(kind.equals("complete")||kind.equals("partial"),"Unknown truth completeness class");validComp(comp);
			require(l.start>0&&l.stop>=l.start&&(!kind.equals("complete")||l.stop-l.start+1>=916),"Complete RF01960 truth requires at least half of CLEN1831");
			final double score=Double.parseDouble(p.parseString(7));require(Double.isFinite(score),"The saved CM reference score must be finite");
			final String key=key(l);require(metadata.get(key)==null,"Duplicate physical CM locus");metadata.put(key,new Meta(kind,comp,score));
			final String old=truthCompartments.get(l.seqid);require(old==null || old.equals(comp),"One truth contig cannot have conflicting compartment labels");truthCompartments.put(l.seqid,comp);
			grade.bucket(l.seqid).truth.add(l);grade.referenceCount++;
		}}finally{require(!in.close(),"Truth read failed");}
		require(!grade.negative||grade.referenceCount==0,"Prok control has zero euk18S truth; do not silently treat another family as18S");
	}
	public void loadCompartments(String path){
		loadCompartments(path,null,-1);
	}
	public void loadCompartments(String path,String catalog,int callerMinLen){
		final ByteFile in=RrnaSourceTargetPrep.open(path);final LineParser1 p=new LineParser1('\t');
		try{header(in,NcrnaStageDiagLogSink.HEADER,path);for(byte[] row;(row=in.nextLine())!=null;){p.set(row);require(p.terms()==21,"Stage sink schema changed");
			if(!p.termEquals("SEED",0)||!p.termEquals("euk18S",1)){continue;}
			final String full=p.parseString(2),seq=RrnaCallerModelAssignment.first(full),comp=compartment(full),old=contigs.get(seq);
			require(old==null||old.equals(comp),"Conflicting compartment for a contig");contigs.put(seq,comp);
		}}finally{require(!in.close(),"Stage header stream read failed");}
		if(catalog!=null){loadSkippedContigs(catalog,callerMinLen);}
		for(RrnaSeedBenchmarkGrade.Bucket b:grade.ordered){
			final String expected=truthCompartments.get(b.seq),observed=contigs.get(b.seq);
			require(observed!=null,"Missing caller diagnostic without bound short-contig evidence: "+b.seq);
			if(!observed.equals(expected)){compartmentAudit.add(new String[]{b.seq,expected,observed,"LABEL_DISAGREEMENT"});}
			// G11 2026-09-29: retained truth labels govern grading; discrepancies
			// are evidence rows, never a reason to abort otherwise valid calls.
			contigs.put(b.seq,expected);
		}
	}
	/** NcrnaScavenger.scavenge returns before emitting SEED for contigs
	 * shorter than minLen. Bind these omissions to sealed extraction metadata;
	 * never invent zero seeds or accept unexplained missing eligible contigs. */
	void loadSkippedContigs(String path,int minimum){
		require(minimum>0,"Skipped-contig evidence requires the actual caller minimum length");
		final ObjectMap<String,RrnaSeedBenchmarkGrade.Bucket> missing=new ObjectMap<String,RrnaSeedBenchmarkGrade.Bucket>(String.class,RrnaSeedBenchmarkGrade.Bucket.class);
		for(RrnaSeedBenchmarkGrade.Bucket b:grade.ordered){if(!b.truth.isEmpty() && contigs.get(b.seq)==null){missing.put(b.seq,b);}}
		if(missing.isEmpty()){return;}
		final ByteFile in=RrnaSourceTargetPrep.open(path);final LineParser1 p=new LineParser1('\t');
		try{header(in,"assembly\tseqid\ttaxid\tlength\tcompartment\tkeyword_class\tfull_header",path);
			for(byte[] row;(row=in.nextLine())!=null;){p.set(row);require(p.terms()==7,"Sealed extraction contig schema changed");if(!p.termEquals(assembly,0)){continue;}
				final String seq=p.parseString(1);final RrnaSeedBenchmarkGrade.Bucket b=missing.get(seq);if(b==null){continue;}
				final int length=p.parseInt(3);final String description=p.parseString(6),comp=compartment(description);validComp(p.parseString(4));
				require(length>0 && length<minimum && RrnaCallerModelAssignment.first(description).equals(seq),"Absent seed diagnostic must bind a genuinely short contig and its actual defline: "+seq);
				require(skippedContigs.add(seq),"Duplicate sealed metadata for skipped contig: "+seq);
				for(Locus locus:b.truth){require(locus.stop<=length,"Skipped-contig truth must lie within the sealed contig length");}
				contigs.put(seq,comp);
				compartmentAudit.add(new String[]{seq,truthCompartments.get(seq),comp,"SHORT_CONTIG_SKIPPED"});
			}
		}finally{require(!in.close(),"Sealed contig metadata read failed");}
		require(skippedContigs.size()==missing.size(),"Every missing diagnostic contig needs independently bound short-contig evidence");
	}
	void writeCompartmentAudit(String prefix){
		RrnaCmWindowParity.fresh(prefix,".compartments.tsv");final ByteStreamWriter out=RrnaSourceTargetPrep.writer(prefix+".compartments.tsv");final ByteBuilder b=new ByteBuilder();
		try{out.println("assembly\tseqid\ttruth_compartment\tobserved_defline_compartment\treason");for(String[] row:compartmentAudit){b.clear().append(assembly);for(String field:row){b.tab().append(field);}out.print(b.nl());}}
		finally{RrnaCmWindowParity.finish(out);}
	}
	public void loadCalls(String path){
		loadCalls(path,true);
	}
	void loadCalls(String path,boolean diagnostics){
		final ByteFile in=RrnaSourceTargetPrep.open(path);final LineParser1 p=new LineParser1('\t');
		try{for(byte[] row;(row=in.nextLine())!=null;){if(row.length==0||row[0]=='#'){continue;}p.set(row);require(p.terms()==9,"Expected GFF9");
			if(!p.termEquals("rRNA",2)&&!p.termEquals("RNA",2)){continue;}
			final String family=RrnaSeedBenchmarkGrade.family(row,p,false);grade.rnaRows++;
			if(family==null){
				final String attributes=p.parseString(8);
				require(!attributes.contains("model:")&&p.termEquals("rRNA",2)&&(attributes.startsWith("5S,")||attributes.startsWith("16S,")||attributes.startsWith("18S,")||attributes.startsWith("23S,")||attributes.startsWith("28S,")),"Unattributed RNA is not a recognized legacy row");
			}else{boolean known=false;for(String name:FAMILIES){if(family.equals(name)){known=true;break;}}require(known,"Unknown RNA family: "+family);}
			if(!"euk18S".equals(family)){grade.otherFamilyRows++;continue;}
			final Locus l=new Locus(p.parseString(0),"euk18S",p.parseInt(3),p.parseInt(4),RrnaSeedBenchmarkGrade.strand(p.parseString(6)));
				require(l.start>0 && l.stop>=l.start,"Caller coordinates must be valid before reciprocal-overlap matching");
				if(diagnostics){require(contigs.get(l.seqid)!=null && !skippedContigs.contains(l.seqid),"Every emitted euk18S call must have an observed seed-stream contig");}
				else{final String comp=truthCompartments.get(l.seqid);contigs.put(l.seqid,comp==null?"unobserved":comp);}
				grade.bucket(l.seqid).calls.add(l);grade.callCount++;
		}}finally{require(!in.close(),"Call GFF read failed");}
		require(grade.callCount+grade.otherFamilyRows==grade.rnaRows,"RNA partitions must conserve output");
	}
	void writeStrata(String prefix,String arm,String clade){
		writeStrata(prefix,arm,clade,COMPS);
	}
	void writeStrata(String prefix,String arm,String clade,String[] compartments){
		final ObjectMap<String,Counts> map=new ObjectMap<String,Counts>(String.class,Counts.class);final ArrayList<Counts> ordered=new ArrayList<Counts>();
		for(String comp:compartments){for(String kind:KINDS){final Counts c=new Counts(comp,kind);map.put(comp+'\t'+kind,c);ordered.add(c);}}
		for(RrnaSeedBenchmarkGrade.Bucket b:grade.ordered){
			for(int i=0;i<b.truth.size();i++){
				final Locus t=b.truth.get(i);final Meta meta=metadata.get(key(t));final int match=b.match.truthMatchedCallIdx[i];
				for(String comp:new String[]{"all",meta.comp}){for(String kind:new String[]{"all",meta.kind}){final Counts c=map.get(comp+'\t'+kind);c.truth++;
					if(match<0){c.fn++;}else{c.tp++;final MatchedPair pair=new MatchedPair("euk18S",t,b.calls.get(match));if(pair.d5==0){c.exact5++;}if(pair.d3==0){c.exact3++;}if(pair.d5==0&&pair.d3==0){c.both++;}}
				}}
			}
			for(int i=0;i<b.calls.size();i++){if(b.match.callMatchedTruthIdx[i]>=0){continue;}final String comp=contigs.get(b.calls.get(i).seqid);map.get("all\tall").fp++;map.get(comp+"\tall").fp++;}
		}
		final Counts all=map.get("all\tall");require(all.tp==grade.tp&&all.fp==grade.fp&&all.fn==grade.fn,"Stratification must preserve the single global one-to-one assignment");
		final ByteStreamWriter out=RrnaSourceTargetPrep.writer(prefix+".strata.tsv");final ByteBuilder b=new ByteBuilder();
		try{out.println(STRATA_HEADER);for(Counts c:ordered){
			b.clear().append(arm).tab().append(clade).tab().append(assembly).tab().append(c.comp).tab().append(c.kind).tab().append(c.truth).tab().append(c.tp).tab();
			if(c.kind.equals("all")){b.append(c.fp);}else{b.append("NA");}
			b.tab().append(c.fn).tab().append(c.exact5).tab().append(c.exact3).tab().append(c.both).tab().append(c.truth>0?1:0).tab().append(c.tp>0?1:0).nl();out.print(b);
		}}finally{RrnaCmWindowParity.finish(out);}
	}
	void writeFailures(String prefix){
		final ByteFile in=RrnaSourceTargetPrep.open(prefix+".stage.loci.tsv.gz");final LineParser1 p=new LineParser1('\t');
		final ByteStreamWriter out=RrnaSourceTargetPrep.writer(prefix+".failures.tsv.gz");final ByteBuilder b=new ByteBuilder();long n=0;
		try{require(in.nextLine()!=null,"Stage reducer must emit header even for zero FNs");out.println("assembly\tseqid\tstart\tstop\tstrand\tcompartment\ttruthClass\tstageEvidence");
			for(byte[] row;(row=in.nextLine())!=null;){p.set(row);require(p.terms()==18&&p.termEquals(assembly,0),"Stage locus binding changed");
				final Locus l=new Locus(p.parseString(1),"euk18S",p.parseInt(2),p.parseInt(3),RrnaSeedBenchmarkGrade.strand(p.parseString(4)));final Meta m=metadata.get(key(l));require(m!=null,"FN must bind one included CM locus");
				String evidence=p.parseString(17);if(skippedContigs.contains(l.seqid)){require(evidence.equals("NO_SEED_OBSERVATION"),"Skipped contig cannot have observed seed/window/path evidence");evidence="CONTIG_BELOW_CALLER_MINLEN";}
				out.print(b.clear().append(assembly).tab().append(l.seqid).tab().append(l.start).tab().append(l.stop).tab().append(l.strand).tab().append(m.comp).tab().append(m.kind).tab().append(evidence).nl());n++;
			}require(n==grade.fn,"Failure strata must conserve all final FNs");
		}finally{require(!in.close(),"Stage loci read failed");RrnaCmWindowParity.finish(out);}
	}
	public static String compartment(String header){require(header!=null,"Compartment needs actual FASTA description");final String h=header.toLowerCase(Locale.ROOT);final boolean m=h.contains("mitochondri")||h.contains("kinetoplast"),p=h.contains("chloroplast")||h.contains("apicoplast")||h.contains("plastid");return m?(p?"ambiguous_both":"mitochondrial"):(p?"plastid":"no_keyword");}
	static void validComp(String c){for(int i=1;i<COMPS.length;i++){if(COMPS[i].equals(c)){return;}}throw new IllegalArgumentException("Unknown compartment: "+c);}
	static String key(Locus l){assert(l!=null):"Physical identity includes contig, both bounds and strand";return l.seqid+'\t'+l.start+'\t'+l.stop+'\t'+l.strand;}
	static void header(ByteFile in,String expected,String path){final byte[] row=in.nextLine();require(row!=null&&new String(row,StandardCharsets.UTF_8).equals(expected),"Unexpected header: "+path);}
	public static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	static final class Meta{Meta(String k,String c,double s){kind=k;comp=c;score=s;}final String kind,comp;final double score;}
	static final class Counts{Counts(String c,String k){comp=c;kind=k;}final String comp,kind;long truth,tp,fp,fn,exact5,exact3,both;}
	final String assembly;public final RrnaSeedBenchmarkGrade grade;
	final ObjectMap<String,Meta> metadata=new ObjectMap<String,Meta>(String.class,Meta.class);
	public final ObjectMap<String,String> contigs=new ObjectMap<String,String>(String.class,String.class);
	final ObjectSet<String> skippedContigs=new ObjectSet<String>(String.class);
	final ObjectMap<String,String> truthCompartments=new ObjectMap<String,String>(String.class,String.class);
	final ArrayList<String[]> compartmentAudit=new ArrayList<String[]>();
	public static final String TRUTH_HEADER="assembly\tseqid\tstart\tstop\tstrand\ttruthClass\tcompartment\tcmScore";
	static final String STRATA_HEADER="arm\tclade\tassembly\tcompartment\ttruthClass\ttruth\tTP\tFP\tFN\texact5\texact3\tbothExact\ttruthBearingGenomes\tcapturedGenomes";
	static final String[] COMPS={"all","no_keyword","mitochondrial","plastid","ambiguous_both"},KINDS={"all","complete","partial"};
	static final String[] WHOLE_COMPS={"all","no_keyword","mitochondrial","plastid","ambiguous_both","unobserved"};
	static final String[] FAMILIES={"euk18S","euk5S","16S","23S","18S","prok5S","r58","lsu","s18","rnasep","srp_small","srp_large","tmrna","sixs_rf00013","sixs_rf01685"};
}
