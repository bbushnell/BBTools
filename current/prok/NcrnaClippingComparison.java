package prok;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.PreParser;
import prok.NcrnaGffFamilyGrader.Locus;
import prok.NcrnaGffFamilyGrader.MatchResult;
import prok.NcrnaGffFamilyGrader.MatchedPair;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/** Saved-call comparison. P1 is inclusive source-annotation-relative; P2/P3
 * use either provisional raw CM hits or reviewed frozen loci with historical masks.
 * The committed producer verifies input and caller hashes before invocation.
 * Reuses maximum-cardinality reciprocal-half matching independently per fold,
 * original assembly/contig/strand. SRP union is a separate cross-check view.
 * @author Ganyu, Raiden */
public final class NcrnaClippingComparison {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args,NcrnaClippingComparison.class,false);
		String refs=null,inputs=null,contigs=null,overlay=null,models=null,runs=null,masks=null,out=null,membership=null,membershipArms=null,truth=null,history=null;int expected=-1;
		long expectedLoci=-1,expectedHistory=-1;int truthVersion=1;
		String baselineSummary=null;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0,"Expected flag=value");
			final String k=arg.substring(0,eq).toLowerCase(Locale.ROOT),v=arg.substring(eq+1);
			if(k.equals("refs")){refs=v;}else if(k.equals("inputs")){inputs=v;}else if(k.equals("contigs")){contigs=v;}
			else if(k.equals("overlay")){overlay=v;}else if(k.equals("models")){models=v;}else if(k.equals("runs")){runs=v;}
			else if(k.equals("masks")){masks=v;}else if(k.equals("outprefix")){out=v;}else if(k.equals("expectedruns")){expected=Integer.parseInt(v);}
			else if(k.equals("membership")){membership=v;}else if(k.equals("membershiparms")){membershipArms=v;}
			else if(k.equals("truth")){truth=v;}else if(k.equals("history")){history=v;}
			else if(k.equals("truthversion")){truthVersion=Integer.parseInt(v);}
			else if(k.equals("baselinesummary")){baselineSummary=v;}
			else if(k.equals("expectedloci")){expectedLoci=Long.parseLong(v);}else if(k.equals("expectedhistoryrows")){expectedHistory=Long.parseLong(v);}
			else{throw new IllegalArgumentException(k);}
		}
		require(refs!=null && inputs!=null && contigs!=null && (overlay!=null ^ truth!=null) && models!=null && runs!=null && masks!=null && out!=null && expected>0,"Explicit inputs and exactly one raw overlay or frozen truth required");
		require(truth==null?history==null:history!=null && expectedLoci>0 && expectedHistory>=0,"Frozen truth requires its historical-mask companion and exact denominators");
		require(truthVersion==1 || truthVersion==2 && truth!=null,"Explicit truthversion=2 requires a frozen version2 reference");
		require(baselineSummary==null || truthVersion==2,"Reference-delta output requires explicit version2 grading");
		Shared.setThreads(1);ReadWrite.setZipThreads(1);
		require((membership==null)==(membershipArms==null),"Explicit subset requires both membership and membershiparms; omitted means the unchanged full-panel contract");
		final NcrnaClippingComparison x=new NcrnaClippingComparison();x.loadModels(models);x.loadReferences(refs);x.loadEukInputs(inputs);x.loadContigs(contigs);
		if(truth==null){x.loadOverlay(overlay);}else{x.loadFrozenTruth(truth,expectedLoci,truthVersion==1?EukPanelTruthDisposition.LEGACY_POLICY:EukPanelTruthFreeze.V2_POLICY,truthVersion);}
		x.loadMasks(masks);if(truth!=null){x.loadHistoricalMasks(history,expectedHistory);}
		if(membership!=null){x.loadMembership(membership,membershipArms);}x.run(runs,out,expected);
		if(baselineSummary!=null){x.writeReferenceDelta(baselineSummary,out);}
		Shared.closeStream(pp.outstream);
	}
	void loadModels(String path){
		assert(modelFamilies.isEmpty()):"The pinned 25-model roster defines RFAM families";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{RrnaResourceIO.header(in,"bundle_sha80\taccession\tfamily\tsubtype\tclen\tga\tsource_cm\tsource_cm_sha80",path);
			for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==8 && known(p.parseString(2)),"Unknown model fold");require(modelFamilies.put(p.parseString(1),p.parseString(2))==null,"Repeated model accession");}
		}finally{require(!in.close(),"Model read failed");}require(modelFamilies.size()==25,"Expected declared25-model roster");
	}
	void loadReferences(String path){
		assert(genomes.isEmpty()):"Load exact source203 references before broader panel";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');int n=0;
		try{RrnaResourceIO.header(in,"index\tgenome\tfasta\tfasta_sha80\tgff\tgff_sha80",path);
			for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==6 && p.parseInt(0)==++n,"Sequential203 source manifest required");
				final Genome g=new Genome(p.parseString(1),"P1",p.parseString(2),p.parseString(3));g.source=p.parseString(4);g.sourceHash=p.parseString(5);addGenome(g);
				RrnaResourceIO.stream(g.fasta,r->{final Contig c=new Contig(g,RrnaResourceIO.first(r.id),r.length(),"PROKARYOTIC");require(g.contigs.put(c.id,c)==null,"Duplicate contig in source genome");g.ordered.add(c);});
				loadSource(g);
			}
		}finally{require(!in.close(),"Source manifest read failed");}require(n==203,"Expected203 source genomes");
	}
	void loadEukInputs(String path){
		assert(genomes.size()==203):"Approved198 panel follows source203";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');int nuclear=0,organelles=0;
		try{p.set(in.nextLine());require(p.terms()==15 && p.termEquals("assembly",1),"Expected bound panel schema");
			for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==15,"Panel schema changed");final boolean nuc=p.termEquals("nuclear",0);
				addGenome(new Genome(p.parseString(1),nuc?"P2":"P3",p.parseString(6),p.parseString(7)));if(nuc){nuclear++;}else{require(p.termEquals("mito",0) || p.termEquals("plastid",0),"Unknown selection group");organelles++;}
			}
		}finally{require(!in.close(),"Panel read failed");}require(nuclear==118 && organelles==80,"Panel scope must remain118+80");
	}
	public void addGenome(Genome g){require(!g.id.isEmpty() && genomes.put(g.id,g)==null,"Assembly identity must be unique across panels: "+g.id);ordered.add(g);}
	void loadContigs(String path){
		assert(genomes.size()==401):"Compartment table binds only the approved198 euk/organelle genomes";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');int n=0;
		try{p.set(in.nextLine());require(p.terms()==15 && p.termEquals("compartment",4),"Expected r6 contig schema");
			for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==15,"Contig schema changed");final Genome g=genomes.get(p.parseString(1));require(g!=null && !g.panel.equals("P1"),"Unknown contig assembly");
				final Contig c=new Contig(g,p.parseString(2),p.parseInt(3),p.parseString(4));require(g.contigs.put(c.id,c)==null,"Duplicate contig identity");g.ordered.add(c);n++;
			}
		}finally{require(!in.close(),"Contig read failed");}require(n==68004,"Expected complete r6 contig inventory");
	}
	public void loadSource(Genome g){
		assert(!g.ordered.isEmpty()):"Actual FASTA lengths precede source coordinate validation";
		final ByteFile in=RrnaResourceIO.open(g.source);final LineParser1 p=new LineParser1('\t');long number=0;
		try{for(byte[] line;(line=in.nextLine())!=null;){number++;if(line.length==0 || line[0]=='#'){continue;}p.set(line);require(p.terms()==9,"Source GFF requires9 fields");
			final boolean region=p.termEquals("region",2);if(!region && !rna(p)){continue;}
			final Contig c=g.contigs.get(p.parseString(0));require(c!=null,"Source contig absent from exact FASTA");
			if(region){final String topology=EukPanelCmProposal.gffAttribute(p.parseString(8),"Is_circular");if(topology!=null){require(topology.equals("true") || topology.equals("false"),"Invalid circular flag");require(p.parseInt(3)==1 && p.parseInt(4)==c.length,"Topology requires full-contig region");c.circular=topology.equals("true");}continue;}
			final String label=sourceLabel(p),family=fold(label);final int start=p.parseInt(3),stop=p.parseInt(4);final String strand=p.parseString(6);
			final Feature f=new Feature(c,family,label,start,stop,strand,g.source+":"+number,p.parseString(8),"SOURCE_ANNOTATION");
			if(!valid(c,start,stop,true) || (!strand.equals("+") && !strand.equals("-"))){f.state="UNSCORABLE_SOURCE";c.unknown.add(f);continue;}
			if(family.equals("UNKNOWN")){f.state="UNKNOWN_SOURCE_FAMILY";c.unknown.add(f);}else{c.truth.add(f);}
		}}
		finally{require(!in.close(),"Source GFF read failed");}
		finishSourceGeometry(g);
	}
	static void finishSourceGeometry(Genome g){
		assert(g!=null):"Resolve source intervals only after all explicit topology rows have been read";
		for(Contig c:g.ordered){
			for(int i=c.truth.size()-1;i>=0;i--){final Feature f=c.truth.get(i);
				if(f.stop>c.length && !c.circular){f.state="UNSCORABLE_SOURCE";c.unknown.add(f);c.truth.remove(i);}
			}
			for(Feature f:c.unknown){if(!valid(c,f.start,f.stop,true) || f.stop>c.length && !c.circular){f.state="UNSCORABLE_SOURCE";}}
			for(Feature f:c.truth){f.geometry=sourceGeometry(c,f);}
			for(Feature f:c.unknown){f.geometry=f.state.equals("UNSCORABLE_SOURCE")?"SOURCE_UNSCORABLE":sourceGeometry(c,f);}
		}
	}
	static String sourceGeometry(Contig c,Feature f){
		assert(valid(c,f.start,f.stop,true) && (f.stop<=c.length || c.circular)):"Physical endpoint labels require a valid source interval";
		if(f.stop>c.length){return "SOURCE_ORIGIN_SPANNING";}
		return f.start==1 || f.stop==c.length?"SOURCE_CONTIG_EDGE":"SOURCE_INTERNAL";
	}
	public String sourceLabel(LineParser1 p){
		assert(p.terms()==9):"Source family attribution uses original GFF fields";
		final String type=p.parseString(2),attrs=p.parseString(8);
		String family=null;boolean unknown=false;final String db=EukPanelCmProposal.gffAttribute(attrs,"Dbxref");
		if(db!=null){for(String token:db.split(",")){if(token.startsWith("RFAM:")){final String f=modelFamilies.get(token.substring(5));if(f==null){unknown=true;}else{if(family!=null && !family.equals(f)){return "UNKNOWN";}family=f;}}}}
		if(unknown){return "UNKNOWN";}
		if(type.equals("tRNA")){return family==null || family.equals("tRNA")?"tRNA":"UNKNOWN";}
		if(type.equals("tmRNA")){return family==null || family.equals("tmRNA")?"tmRNA":"UNKNOWN";}
		if(family!=null){return family;}
		if(type.equals("rRNA")){if(attrs.contains("5.8S")){return "5.8S";}final String old=LegacyRrnaGffGrader.extractLegacyFamily(attrs,p.parseInt(4)-p.parseInt(3)+1);if(old!=null){return old;}}
		return "UNKNOWN";
	}
	void loadOverlay(String path){
		assert(rawHits==0):"Retain every raw hit ID; do not collapse model hits into biological truth";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');final ObjectSet<String> seen=new ObjectSet<String>(String.class);
		try{RrnaResourceIO.header(in,EukPanelLoci.HEADER,path);for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==40,"Expected40-column raw review overlay");
			final Genome g=genomes.get(p.parseString(0));require(g!=null && !g.panel.equals("P1"),"Raw hit outside198 panel");final Contig c=g.contigs.get(p.parseString(10));
			require(c!=null && p.termEquals(c.compartment,11) && p.termEquals("VALID_INTERVAL",13),"Hit must retain actual valid contig/compartment");
			final String family=p.parseString(7),id=p.parseString(1)+":"+p.parseLong(4);require(known(family) && family.equals(modelFamilies.get(p.parseString(6))) && seen.add(id),"Unknown/duplicate raw hit");
			final int from=p.parseInt(19),to=p.parseInt(20);final Feature f=new Feature(c,family,p.parseString(6),Math.min(from,to),Math.max(from,to),p.parseString(21),id,"NA","RAW_GA_HITS_PROVISIONAL");
			f.truncated=!p.termEquals("no",22);
			require(valid(c,f.start,f.stop,false),"Raw hit outside contig");f.geometry=!p.termEquals("INTERNAL",31) || !p.termEquals("no",22)?"EDGE_OR_TRUNCATED":p.termEquals("true",32)?"FULLSPAN_INTERNAL":"PARTIAL_INTERNAL";
			if(p.parseDouble(24)>=p.parseDouble(26)){require(p.termEquals("true",35),"GA flag disagrees");c.truth.add(f);}else{require(p.termEquals("false",35),"GA flag disagrees");f.state="BELOW_GA";c.masks.add(f);}rawHits++;
		}}finally{require(!in.close(),"Overlay read failed");}require(rawHits==88223,"Raw-hit denominator changed");
	}
	/** Frozen mode retains the original matching algorithm; only the reference unit
	 * and declared metric eligibility change. Raw-overlay mode stays byte-compatible. */
	public void loadFrozenTruth(String path,long expected){
		loadFrozenTruth(path,expected,EukPanelTruthDisposition.LEGACY_POLICY);
	}
	public void loadFrozenTruth(String path,long expected,String policy){
		loadFrozenTruth(path,expected,policy,1);
	}
	public void loadFrozenTruth(String path,long expected,String policy,int version){
		assert(rawHits==0 && referenceRows==0):"Raw hits and frozen loci are mutually exclusive reference universes";
		EukPanelTruthFreeze.validateVersion(policy,version);
		final boolean prok=policy.equals(EukPanelTruthDisposition.PROK_POLICY);
		require(expected>0 && !frozenMode,"One complete frozen truth table required");frozenMode=true;expectedFrozenLoci=expected;
		frozenReferenceKind="FROZEN_LOCUS_TRUTH_V"+version;
		for(Genome g:ordered){if(!g.panel.equals("P1")){g.referenceKind=prok?"FRESH_PROK_CM_TRUTH_V1":frozenReferenceKind;}}
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');final ObjectSet<String> seen=new ObjectSet<String>(String.class);
		try{RrnaResourceIO.header(in,EukPanelTruthFreeze.HEADER,path);for(byte[] line;(line=in.nextLine())!=null;){p.set(line);
			require(p.terms()==35 && p.termEquals("TRUTH_V"+version,32) && p.termEquals(policy,34),"Only the explicitly selected frozen truth policy/version is supported");
			final Genome g=genomes.get(p.parseString(1));require(g!=null && !g.panel.equals("P1"),"Frozen euk truth must bind a selected euk/organelle genome");
			final Contig c=g.contigs.get(p.parseString(2));require(c!=null && p.termEquals(c.compartment,8),"Frozen contig/compartment binding changed");
			if(prok){require(c.compartment.equals("PROKARYOTIC") && (g.panel.equals("bacteria") || g.panel.equals("archaea")) && p.termEquals(g.panel,7),"Fresh prok truth must retain its domain and physical compartment");}
			final String id=p.parseString(0),family=p.parseString(6),seed=p.parseString(9),subtype=p.parseString(25),admission=p.parseString(26),flags=p.parseString(28);
			require(seen.add(id) && known(family) && family.equals(modelFamilies.get(seed)),"Duplicate locus or unknown seed/fold binding");
			final boolean ga=EukPanelLocusDiagnostics.bit(p,14),seedGa=EukPanelLocusDiagnostics.bit(p,24);
			require(admission.equals(ga?"ADMITTED":"BELOW_GA"),"Frozen admission must retain any-GA policy; withdrawn masks cannot appear");
			require(ga?family.equals(modelFamilies.get(subtype)):subtype.equals("NA"),"Frozen subtype must be a same-fold supported member, or NA below GA");
			final Feature f=new Feature(c,family,ga?subtype:seed,p.parseInt(4),p.parseInt(5),p.parseString(3),"LOCUS:"+id,flags,frozenReferenceKind);
			require(valid(c,f.start,f.stop,false) && (f.strand.equals("+") || f.strand.equals("-")),"Frozen locus outside its physical contig");
			f.geometry=p.parseString(27);require(f.geometry.equals("FULLSPAN_INTACT") || f.geometry.equals("EDGE") || f.geometry.equals("PARTIAL_INTERNAL") || f.geometry.equals("UNASSESSABLE"),"Unknown frozen geometry stratum");
			f.truncated=EukPanelLocusDiagnostics.bit(p,16);f.boundaryEligible=EukPanelLocusDiagnostics.bit(p,31);f.bothEligible=EukPanelLocusDiagnostics.bit(p,33);
			require(f.boundaryEligible==(ga && seedGa) && f.bothEligible==(f.boundaryEligible && !EukPanel58sHomologAudit.member(flags,"AMBIGUOUS_CROSSFOLD")),"Frozen boundary eligibility must preserve seed support and both-only crossfold exclusion");
			if(ga){c.truth.add(f);}else{f.state="BELOW_GA";c.masks.add(f);}referenceRows++;
		}}finally{require(!in.close(),"Frozen truth read failed");}require(referenceRows==expected,"Frozen locus denominator changed");
	}
	/** History is a separate retained source inventory, never new TP/FN truth.
	 * Eligible same-fold masks run only after admitted matching and duplicates. */
	public void loadHistoricalMasks(String path,long expected){
		assert(frozenMode && referenceRows>0):"Historical masks supplement an explicitly frozen reference set";
		require(!historyLoaded && expected>=0,"One provenance-bound historical inventory required");
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');final ObjectSet<String> seen=new ObjectSet<String>(String.class);
		try{RrnaResourceIO.header(in,HISTORY_HEADER,path);for(byte[] line;(line=in.nextLine())!=null;){p.set(line);
			require(p.terms()==14,"Historical-mask schema changed");final String assembly=p.parseString(0),source=p.parseString(8);final long sourceLine=p.parseLong(10);
			require(sourceLine>0 && SeedOffsetTable.validSha80(p.parseString(9)) && seen.add(assembly+'\t'+source+'\t'+sourceLine),"Historical source identities must be unique and provenance-bound");
			final Genome g=genomes.get(assembly);require(g!=null && !g.panel.equals("P1"),"Historical euk masks cannot replace prok source truth");historyRows++;
			final boolean eligible=EukPanelLocusDiagnostics.bit(p,7);if(!eligible){continue;}
			final String state=p.parseString(6),family=p.parseString(5);require(state.equals("HIST_ONLY") || state.equals("HIST_LINKED_UNSCORED"),"Only unresolved historical RNA may mask an unmatched call");
			require(known(family),"Unclassified historical RNA cannot mask arbitrary folds");final Contig c=g.contigs.get(p.parseString(1));require(c!=null,"Eligible historical mask has unresolved contig");
			final Feature f=new Feature(c,family,family,p.parseInt(3),p.parseInt(4),p.parseString(2),source+":"+sourceLine,p.parseString(12),"HIST_ONLY");
			require(valid(c,f.start,f.stop,false) && (f.strand.equals("+") || f.strand.equals("-")),"Eligible historical mask must be an intact same-frame source interval");
			f.state=state;c.masks.add(f);historicalMasks++;
		}}finally{require(!in.close(),"Historical-mask read failed");}require(historyRows==expected,"Historical inventory denominator changed");historyLoaded=true;
	}
	void loadMasks(String path){
		assert(rawHits>0 || frozenMode && referenceRows>0):"Outside-map mask supplements loaded reference evidence, never changes it";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{byte[] line=in.nextLine();require(line!=null && line[0]=='#',"Mask policy header required");p.set(in.nextLine());require(p.terms()==33 && p.termEquals("mask_eligible",31),"Expected reviewed source-mask schema");
			for(;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==33 && p.termEquals("SRP",29),"Only explicitly union-scoped SRP masks supported");if(!p.termEquals("1",31)){continue;}
				final Genome g=genomes.get(p.parseString(1));require(g!=null,"Mask assembly missing");final Contig c=g.contigs.get(p.parseString(8));require(c!=null,"Mask contig missing");
				final Feature f=new Feature(c,"SRP",p.parseString(30),p.parseInt(17),p.parseInt(18),p.parseString(20),p.parseString(2)+":"+p.parseLong(4),p.parseString(22),"OUT_OF_MAP");require(valid(c,f.start,f.stop,false),"Mask outside contig");f.state="OUT_OF_MAP";c.masks.add(f);
			}
		}finally{require(!in.close(),"Mask read failed");}
	}
	void run(String path,String prefix,int expected){
		assert(frozenMode?referenceRows==expectedFrozenLoci && historyLoaded:rawHits==88223):"Grade only after complete reference and required historical-mask loading";
		RrnaResourceIO.fresh(prefix,".summary.tsv",".truth_strata.tsv",".emitted_labels.tsv",".audit.tsv.gz",".unknown.tsv.gz",".runtime.tsv",".call_overlaps.tsv.gz",".workload.tsv",".physical.tsv");
		final ByteStreamWriter summary=RrnaResourceIO.writer(prefix+".summary.tsv"),strata=RrnaResourceIO.writer(prefix+".truth_strata.tsv"),labels=RrnaResourceIO.writer(prefix+".emitted_labels.tsv");
		final ByteStreamWriter audit=RrnaResourceIO.writer(prefix+".audit.tsv.gz"),unknown=RrnaResourceIO.writer(prefix+".unknown.tsv.gz"),runtime=RrnaResourceIO.writer(prefix+".runtime.tsv");
		final ByteStreamWriter overlaps=RrnaResourceIO.writer(prefix+".call_overlaps.tsv.gz");
		final ByteStreamWriter workload=RrnaResourceIO.writer(prefix+".workload.tsv");
		final ByteStreamWriter physical=RrnaResourceIO.writer(prefix+".physical.tsv");
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');final ByteBuilder b=new ByteBuilder();final ObjectSet<String> seen=new ObjectSet<String>(String.class);int n=0;
		try{
			audit.println("arm\tpanel\tassembly\tcontig\tcompartment\tview\trole\tid\tlabel\tstart\tstop\tstrand\toutcome\tpartner\td5\td3\tgeometry\treference_kind");
			unknown.println("arm\tpanel\tassembly\tcontig\trole\tid\tstart\tstop\tstrand\treason\toriginal_attributes");
			runtime.println("arm\tpanel\ttask\tgenomes\tprocess_wall_seconds\tuser_seconds\tsystem_seconds");
			workload.println("# Caller counters are reported as logged, not final GFF counts; absent families are not inferred as zero. GLOBAL metrics must not be summed with family alignments.");
			workload.println("arm\tpanel\ttask\tmetric\tfamily\tcount");
			overlaps.println("arm\tpanel\tassembly\tcontig\tview\tcall_a\tcall_b\toverlap_bp\texact_coordinates\tstatus");
			for(Genome g:ordered){if(membershipPath!=null && !membershipMembers.contains(g.panel+'\t'+g.id)){continue;}for(Contig c:g.ordered){for(Feature f:c.unknown){unknown(unknown,b,"REFERENCE",f,"SOURCE",f.state);}}}
			RrnaResourceIO.header(in,"arm\tpanel\ttask\tresult_root\tseal_sha80\tcalls_gff\tcalls_sha80\tinput_genomes",path);
			for(byte[] line;(line=in.nextLine())!=null;){
				p.set(line);require(p.terms()==8,"Verified run schema changed");final String arm=p.parseString(0),panel=p.parseString(1),task=p.parseString(2),root=p.parseString(3),calls=p.parseString(5);
				require(seen.add(arm+'\t'+panel+'\t'+task) && SeedOffsetTable.validSha80(p.parseString(4)) && SeedOffsetTable.validSha80(p.parseString(6)),"Duplicate/unsealed caller run");
				final ArrayList<Genome> selected=runInputs(root+"/inputs.tsv",panel,p.parseInt(7));final ObjectMap<String,Contig> byContig=new ObjectMap<String,Contig>(String.class,Contig.class);
				registerCoverage(arm,panel,selected);
				for(Genome g:selected){for(Contig c:g.ordered){require(c.calls.isEmpty() && byContig.put(c.id,c)==null,"Combined run has ambiguous repeated contig IDs");}}
				loadCalls(calls,arm,byContig,unknown,b);
				for(Genome g:selected){for(Contig c:g.ordered){
					for(String view:FAMILIES){grade(c,arm,view,audit,b);callOverlaps(c,arm,view,overlaps,b);}
					c.calls.clear();
				}}
				readTime(root+"/time.txt",runtime,b,arm,panel,task,selected.size());
				readWorkload(root+"/run.stderr",workload,b,arm,panel,task);n++;
			}
			require(n==expected,"Verified caller-run denominator changed: "+n+" != "+expected);
			validateCoverage();
			if(membershipPath!=null){summary.println("# EXPLICIT_SUBSET: all declared arms cover identical membership; not a full-panel estimate. See scope.tsv sidecar; membership="+membershipPath);}
			if(frozenMode){summary.println("# P1=SOURCE_ANNOTATION_RELATIVE inclusive; P2/P3="+frozenReferenceKind+". HIST_ONLY masks applied after one-to-one matching and duplicate handling. Precision/FP are pooled by fold and compartment; geometry rows have reference recall and boundary metrics, not assigned FPs. Crossfold excludes only both-end metrics; boundary fractions use explicit eligible-TP denominators. PARTIAL_INTERNAL remains included in labelled pooled rows and separate in strata.");}
			else{summary.println("# P1=SOURCE_ANNOTATION_RELATIVE inclusive; P2/P3=RAW_GA_HITS_PROVISIONAL, NOT unique genes; HIST_ONLY_MASK_NOT_APPLIED. Below-GA/outside-map masks applied. SRP union is a separate overlapping view. Apparent FP/FN counts are reference-relative, not confirmed biological errors. Raw duplicate/overlapping predictions are listed independently of reference matching.");}
			summary.println("arm\tpanel\tcompartment\tview\tTP\tapparentFP\tFN\texact5\texact3\tbothExact\tduplicateUnmatchedCalls\tunscoredCalls"+(frozenMode?"\tboundaryEligibleTP\tbothEndsEligibleTP":"")+"\tprecision\trecall\tbothExactFraction"+(frozenMode?"\texact5Fraction\texact3Fraction":""));
			for(String key:EukPanelCmProposal.sorted(totals)){final long[] v=totals.get(key);b.clear().append(key);for(long value:v){b.tab().append(value);}fraction(b,v[0],v[0]+v[1]);fraction(b,v[0],v[0]+v[2]);fraction(b,v[5],frozenMode?v[9]:v[0]);if(frozenMode){fraction(b,v[3],v[8]);fraction(b,v[4],v[8]);}summary.print(b.nl());}
			strata.println("arm\tpanel\tcompartment\tview\ttruth_geometry\tTP\tFN\texact5\texact3\tbothExact"+(frozenMode?"\tboundaryEligibleTP\tbothEndsEligibleTP":""));
			for(String key:EukPanelCmProposal.sorted(truthStrata)){b.clear().append(key);for(long value:truthStrata.get(key)){b.tab().append(value);}strata.print(b.nl());}
			labels.println("arm\tpanel\tcompartment\temitted_label\tcalls\tdelta_vs_A\tdelta_vs_D");
			for(String ap:EukPanelCmProposal.sorted(armPanels)){final int tab=ap.indexOf('\t');final String arm=ap.substring(0,tab),panel=ap.substring(tab+1);
				for(String key:EukPanelCmProposal.sorted(labelUniverse)){if(!key.startsWith(panel+'\t')){continue;}final long count=EukPanelCmProposal.get(labelCounts,arm+'\t'+key);
					b.clear().append(arm).tab().append(key).tab().append(count);
					for(String control:new String[]{"A","D"}){b.tab();if(armPanels.get(control+'\t'+panel)==null){b.append("NA");}else{b.append(count-EukPanelCmProposal.get(labelCounts,control+'\t'+key));}}labels.print(b.nl());
				}
			}
			physical.println("arm\tpanel\tcompartment\tview\tbasis\tcontact\treference_truncated\tTP\tapparentFP\tFN\tunscoredCalls");
			for(String key:EukPanelCmProposal.sorted(physicalCounts)){b.clear().append(key);for(long value:physicalCounts.get(key)){b.tab().append(value);}physical.print(b.nl());}
		}finally{final boolean error=in.close();RrnaResourceIO.finish(summary,strata,labels,audit,unknown,runtime,overlaps,workload,physical);require(!error,"Run roster read failed");}
		if(frozenMode){if(membershipPath!=null){writeScope(prefix);}System.err.println("FROZEN_LOCUS_COMPARISON_COMPLETE runs="+n+" reference_genomes=401 reference_loci="+referenceRows+" history_rows="+historyRows+" historical_masks="+historicalMasks+" no_caller_runs=true");}
		else if(membershipPath!=null){writeScope(prefix);System.err.println("CLIPPING_COMPARISON_SUBSET_COMPLETE runs="+n+" reference_genomes=401 evaluated_genomes_per_arm="+membershipMembers.size()+" arms="+membershipArms.size()+" raw_hits="+rawHits+" provisional=true no_truth_freeze=true");}
		else{System.err.println("CLIPPING_COMPARISON_COMPLETE runs="+n+" genomes=401 raw_hits="+rawHits+" provisional=true no_truth_freeze=true");}
	}
	/** Pair native count summaries from the same saved caller arms. No rematching.
	 * P1 source truth is untouched by the eukaryotic addendum and must remain exact. */
	public void writeReferenceDelta(String baseline,String prefix){
		assert(frozenMode && !totals.isEmpty()):"Reference deltas follow complete frozen-reference grading";
		final ObjectMap<String,long[]> before=new ObjectMap<String,long[]>(String.class,long[].class);
		final ByteFile in=RrnaResourceIO.open(baseline);final LineParser1 p=new LineParser1('\t');boolean header=false;
		try{for(byte[] line;(line=in.nextLine())!=null;){if(line.length>0 && line[0]=='#'){continue;}p.set(line);
			if(!header){require(p.terms()==19 && p.termEquals("TP",4) && p.termEquals("boundaryEligibleTP",12),"Expected native frozen-reference19-column summary");header=true;continue;}
			require(p.terms()==19,"Baseline summary schema changed");
			final String key=p.parseString(0)+'\t'+p.parseString(1)+'\t'+p.parseString(2)+'\t'+p.parseString(3);
			final long[] old=new long[10],now=totals.get(key);
			for(int i=0;i<old.length;i++){old[i]=p.parseLong(i+4);require(old[i]>=0,"Baseline metric counts cannot be negative");}
			require(now!=null && before.put(key,old)==null,"Delta requires identical unique arm/panel/compartment/family groups");
			if(p.termEquals("P1",1)){require(java.util.Arrays.equals(old,now),"Eukaryotic reference revision changed the untouched prokaryotic panel: "+key);}
		}}finally{require(!in.close(),"Baseline summary read failed");}
		require(header && before.size()==totals.size(),"Every revised summary group needs a baseline counterpart");
		RrnaResourceIO.fresh(prefix,".reference_delta.tsv");final ByteStreamWriter out=RrnaResourceIO.writer(prefix+".reference_delta.tsv");final ByteBuilder b=new ByteBuilder();
		try{out.println("arm\tpanel\tcompartment\tview\tbaseline_TP\tbaseline_FP\tbaseline_FN\trevised_TP\trevised_FP\trevised_FN\tdelta_TP\tdelta_FP\tdelta_FN\tbaseline_unscored\trevised_unscored\tbaseline_precision\trevised_precision\tbaseline_recall\trevised_recall\tbaseline_both_exact_fraction\trevised_both_exact_fraction");
			for(String key:EukPanelCmProposal.sorted(before)){final long[] old=before.get(key),now=totals.get(key);b.clear().append(key);
				for(int i=0;i<3;i++){b.tab().append(old[i]);}for(int i=0;i<3;i++){b.tab().append(now[i]);}for(int i=0;i<3;i++){b.tab().append(now[i]-old[i]);}
				b.tab().append(old[7]).tab().append(now[7]);fraction(b,old[0],old[0]+old[1]);fraction(b,now[0],now[0]+now[1]);
				fraction(b,old[0],old[0]+old[2]);fraction(b,now[0],now[0]+now[2]);fraction(b,old[5],old[9]);fraction(b,now[5],now[9]);out.print(b.nl());
			}
		}finally{RrnaResourceIO.finish(out);}
	}
	/** Scope is an explicit immutable input, not inferred from whichever runs finish.
	 * The producer pins its digest alongside the caller and reference manifests. */
	public void loadMembership(String path,String arms){
		assert(!genomes.isEmpty() && covered.isEmpty()):"Subset membership must bind existing reference identities before any caller run";
		require(membershipPath==null && path!=null && arms!=null,"Load exactly one explicit subset and expected-arm set");
		for(String arm:arms.split(",",-1)){require(arm.matches("[A-Za-z0-9_]+") && membershipArms.add(arm),"Subset arms must be nonempty unique safe identifiers");}
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{RrnaResourceIO.header(in,"panel\tassembly",path);for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==2,"Subset membership requires panel and assembly only");
			final String panel=p.parseString(0),id=p.parseString(1);final Genome g=genomes.get(id);
			require(g!=null && g.panel.equals(panel),"Subset member must match an existing reference panel/assembly: "+panel+'/'+id);
			require(membershipMembers.add(panel+'\t'+id),"Duplicate subset member: "+panel+'/'+id);EukPanelCmProposal.increment(membershipPanels,panel,1);
		}}finally{require(!in.close(),"Subset membership read failed");}
		require(!membershipMembers.isEmpty(),"Empty subset is not an evaluation scope");membershipPath=path;
	}
	public void writeScope(String prefix){
		assert(membershipPath!=null):"Only explicit-subset output needs a separate membership sidecar";
		RrnaResourceIO.fresh(prefix,".scope.tsv");final ByteStreamWriter out=RrnaResourceIO.writer(prefix+".scope.tsv");final ByteBuilder b=new ByteBuilder();
		try{out.println("scope\tarm\tpanel\tassembly");for(String arm:sorted(membershipArms)){for(String member:sorted(membershipMembers)){out.print(b.clear().append("EXPLICIT_SUBSET").tab().append(arm).tab().append(member).nl());}}}
		finally{RrnaResourceIO.finish(out);}
	}
	static ArrayList<String> sorted(ObjectSet<String> set){
		assert(set!=null):"Scope output must be deterministic independently of hash iteration";
		final ArrayList<String> list=new ArrayList<String>();for(String key:set.keys()){if(key!=null){list.add(key);}}Collections.sort(list);return list;
	}
	public void registerCoverage(String arm,String panel,ArrayList<Genome> selected){
		assert(selected!=null && !selected.isEmpty()):"Every represented arm-panel must consume explicit genomes";
		if(membershipPath!=null){require(membershipArms.contains(arm) && membershipPanels.get(panel)!=null,"Undeclared subset arm/panel: "+arm+'/'+panel);}
		final String ap=arm+'\t'+panel;EukPanelCmProposal.increment(armPanels,ap,1);
		for(Genome g:selected){
			if(membershipPath!=null){require(membershipMembers.contains(panel+'\t'+g.id),"Undeclared subset genome: "+panel+'/'+g.id);}
			require(g.panel.equals(panel) && covered.add(ap+'\t'+g.id),"Repeated arm/panel/genome or wrong panel: "+ap+'/'+g.id);EukPanelCmProposal.increment(coverageCounts,ap,1);}
	}
	public void validateCoverage(){
		assert(!armPanels.isEmpty()):"Completeness is checked per represented arm-panel, not inferred from total run rows";
		if(membershipPath!=null){
			for(String arm:sorted(membershipArms)){
				for(String panel:EukPanelCmProposal.sorted(membershipPanels)){final String ap=arm+'\t'+panel;
					require(armPanels.get(ap)!=null && EukPanelCmProposal.get(coverageCounts,ap)==membershipPanels.get(panel)[0],"Incomplete subset arm/panel coverage: "+ap);}
				for(String member:sorted(membershipMembers)){require(covered.contains(arm+'\t'+member),"Missing expected subset member: "+arm+'/'+member);}
			}
			require(covered.size()==(long)membershipArms.size()*membershipMembers.size(),"Subset coverage must conserve each declared arm/member exactly once");return;
		}
		for(String ap:EukPanelCmProposal.sorted(armPanels)){
			final String panel=ap.substring(ap.indexOf('\t')+1);final int expected=panel.equals("P1")?203:panel.equals("P2")?118:panel.equals("P3")?80:-1;
			require(expected>0 && EukPanelCmProposal.get(coverageCounts,ap)==expected,"Incomplete arm/panel genome coverage: "+ap);
			int found=0;for(Genome g:ordered){if(g.panel.equals(panel)){require(covered.contains(ap+'\t'+g.id),"Missing expected arm/panel/genome: "+ap+'/'+g.id);found++;}}
			require(found==expected,"Reference panel membership changed: "+panel);
		}
	}
	ArrayList<Genome> runInputs(String path,String panel,int expected){
		assert(expected>0):"A caller run consumes an explicit nonempty genome roster";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');final ArrayList<Genome> selected=new ArrayList<Genome>();final ObjectSet<String> seen=new ObjectSet<String>(String.class);
		try{RrnaResourceIO.header(in,"index\tassembly\tfasta\tsha80",path);for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==4,"Caller input schema changed");final Genome g=genomes.get(p.parseString(1));
			require(g!=null && g.panel.equals(panel) && p.termEquals(g.fasta,2) && seen.add(g.id) && SeedOffsetTable.validSha80(p.parseString(3)),"Caller input identity must match its reference panel");
			require(g.hash.equals("NA") || p.termEquals(g.hash,3),"Caller/reference FASTA binding changed");selected.add(g);
		}}finally{require(!in.close(),"Caller input roster read failed");}require(selected.size()==expected,"Caller genome count changed");return selected;
	}
	public void loadCalls(String path,String arm,ObjectMap<String,Contig> byContig,ByteStreamWriter unknown,ByteBuilder b){
		assert(!byContig.isEmpty()):"Caller coordinates must resolve within this run, not another assembly";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');long number=0,rnaRows=0,classified=0,unknownRows=0;
		try{for(byte[] line;(line=in.nextLine())!=null;){number++;if(line.length==0 || line[0]=='#'){continue;}p.set(line);require(p.terms()==9,"Caller GFF requires9 columns");if(!rna(p)){continue;}
			rnaRows++;
			final Contig c=byContig.get(p.parseString(0));require(c!=null,"Caller contig outside declared run");final int start=p.parseInt(3),stop=p.parseInt(4);final String strand=p.parseString(6),attrs=p.parseString(8);
			require(valid(c,start,stop,false) && (strand.equals("+") || strand.equals("-")),"Invalid caller RNA interval/strand");
			String label=PanelAnnotationGrade.colonField(attrs,"family");if(label==null){label=p.termEquals("tRNA",2)?"tRNA":p.termEquals("tmRNA",2)?"tmRNA":LegacyRrnaGffGrader.extractLegacyFamily(attrs,stop-start+1);}
			if(label==null){label="UNKNOWN";}final Feature f=new Feature(c,fold(label),label,start,stop,strand,path+":"+number,attrs,"CALL");
			EukPanelCmProposal.increment(labelCounts,arm+'\t'+c.genome.panel+'\t'+c.compartment+'\t'+label,1);
			EukPanelCmProposal.increment(labelUniverse,c.genome.panel+'\t'+c.compartment+'\t'+label,1);
			if(f.family.equals("UNKNOWN")){unknown(unknown,b,arm,f,"CALL","UNKNOWN_CALL_FAMILY");unknownRows++;}else{c.calls.add(f);classified++;}
		}}finally{require(!in.close(),"Caller GFF read failed");}
		require(rnaRows==classified+unknownRows,"Every caller RNA row must be classified or explicitly diagnosed");
		require(unknownRows==0,"Unknown caller RNA family rows="+unknownRows+" in "+path+"; diagnostic rows retained, refusing incomplete scored denominators");
	}
	public void grade(Contig c,String arm,String view,ByteStreamWriter audit,ByteBuilder b){
		assert(known(view)):"Only declared biological folds and explicit SRP union are graded";
		final ArrayList<Feature> truths=new ArrayList<Feature>(),calls=new ArrayList<Feature>();final ArrayList<Locus> tl=new ArrayList<Locus>(),cl=new ArrayList<Locus>();
		for(Feature f:c.truth){if(inView(f,view)){truths.add(f);tl.add(f.locus(view));}}for(Feature f:c.calls){if(inView(f,view)){calls.add(f);cl.add(f.locus(view));}}
		final int circle=c.circular?c.length:0;final MatchResult match=NcrnaGffFamilyGrader.matchOneToOne(cl,tl,.5f,circle);
		final String key=arm+'\t'+c.genome.panel+'\t'+c.compartment+'\t'+view;long[] v=totals.get(key);if(v==null){v=new long[frozenMode?10:8];totals.put(key,v);}
		for(int i=0;i<truths.size();i++){
			final Feature t=truths.get(i);final int j=match.truthMatchedCallIdx[i];final Feature call=j<0?null:calls.get(j);int d5=0,d3=0;
			final String skey=key+'\t'+t.geometry;long[] s=truthStrata.get(skey);if(s==null){s=new long[frozenMode?7:5];truthStrata.put(skey,s);}
			if(call==null){v[2]++;s[1]++;}else{
				final MatchedPair pair=new MatchedPair(view,tl.get(i),cl.get(j));d5=delta(pair.d5,circle);d3=delta(pair.d3,circle);v[0]++;s[0]++;
				if(t.boundaryEligible){if(frozenMode){v[8]++;s[5]++;}if(d5==0){v[3]++;s[2]++;}if(d3==0){v[4]++;s[3]++;}}
				if(t.bothEligible){if(frozenMode){v[9]++;s[6]++;}if(d5==0 && d3==0){v[5]++;s[4]++;}}
			}
			row(audit,b,arm,view,"TRUTH",t,call==null?"FN":"TP",call,d5,d3);
			if(!c.genome.panel.equals("P1")){physical(key,t,"REFERENCE",call==null?2:0);}
		}
		for(int i=0;i<calls.size();i++){
			final Feature call=calls.get(i);final int j=match.callMatchedTruthIdx[i];final Feature truth=j<0?null:truths.get(j);Feature relation=truth;String status="TP";
			if(truth==null){
				if(match.callIsDuplicate[i]){status="DUPLICATE_APPARENT_FP";v[1]++;v[6]++;}
				else{relation=mask(c,cl.get(i),view,circle);if(relation!=null){status=relation.origin.equals("HIST_ONLY")?"UNSCORED_HIST_ONLY":"UNSCORED_OVERLAP";v[7]++;}else{status="APPARENT_FP";v[1]++;}}
			}
			int d5=0,d3=0;if(truth!=null){final MatchedPair pair=new MatchedPair(view,tl.get(j),cl.get(i));d5=delta(pair.d5,circle);d3=delta(pair.d3,circle);}
			row(audit,b,arm,view,"CALL",call,status,relation,d5,d3);
			if(!c.genome.panel.equals("P1") && truth==null){physical(key,call,"CALL",status.startsWith("UNSCORED_")?3:1);}
		}
	}
	void physical(String key,Feature f,String basis,int column){
		assert(column>=0 && column<4):"Physical counters retain TP/FP/FN/unscored as distinct outcomes";
		final String contact=f.start==1 || f.stop==f.contig.length?"EDGE":"INTERIOR";
		final String k=key+'\t'+basis+'\t'+contact+'\t'+(basis.equals("REFERENCE")?(f.truncated?"YES":"NO"):"NA");
		long[] v=physicalCounts.get(k);if(v==null){v=new long[4];physicalCounts.put(k,v);}v[column]++;
	}
	public static void callOverlaps(Contig c,String arm,String view,ByteStreamWriter out,ByteBuilder b){
		assert(known(view)):"Raw-call overlap is a separate diagnostic, never deduplication or truth matching";
		final ArrayList<Feature> sorted=new ArrayList<Feature>();for(Feature f:c.calls){if(inView(f,view)){require(valid(c,f.start,f.stop,false),"Caller-overlap sweep requires non-wrapping valid intervals");sorted.add(f);}}
		Collections.sort(sorted,CALL_ORDER);
		for(int i=0;i<sorted.size();i++){final Feature a=sorted.get(i);
			for(int j=i+1;j<sorted.size();j++){final Feature z=sorted.get(j);if(!a.strand.equals(z.strand) || z.start>a.stop){break;}
				final int n=c.circular?NcrnaGffFamilyGrader.circularOverlap(a.locus(view),z.locus(view),c.length):Math.min(a.stop,z.stop)-Math.max(a.start,z.start)+1;if(n<=0){continue;}
				out.print(b.clear().append(arm).tab().append(c.genome.panel).tab().append(c.genome.id).tab().append(c.id).tab().append(view).tab().append(a.id).tab().append(z.id).tab().append(n)
					.tab().append(a.start==z.start && a.stop==z.stop?1:0).tab().append("RAW_CALL_DIAGNOSTIC_NOT_DEDUPLICATED").nl());
			}
		}
	}
	public static Feature mask(Contig c,Locus call,String view,int circle){
		assert(call.strand=='+' || call.strand=='-'):"Unscored overlap still requires a supported caller strand";
		for(Feature f:c.masks){if(inView(f,view) || f.family.equals("SRP") && view.startsWith("SRP")){if(overlap(call,f,view,circle)){return f;}}}
		for(Feature f:c.unknown){if(f.state.equals("UNKNOWN_SOURCE_FAMILY") && overlap(call,f,view,circle)){return f;}}
		return null;
	}
	static boolean overlap(Locus call,Feature f,String view,int circle){
		if(!f.strand.equals(Character.toString(call.strand))){return false;}
		final Locus truth=f.locus(view);final int n=circle==0?Math.min(call.stop,truth.stop)-Math.max(call.start,truth.start)+1:NcrnaGffFamilyGrader.circularOverlap(call,truth,circle);
		return n>0 && n>=.5*call.length() && n>=.5*truth.length();
	}
	static void row(ByteStreamWriter out,ByteBuilder b,String arm,String view,String role,Feature f,String outcome,Feature partner,int d5,int d3){
		assert(f!=null):"Each audit row identifies one original source or call";
		b.clear().append(arm).tab().append(f.contig.genome.panel).tab().append(f.contig.genome.id).tab().append(f.contig.id).tab().append(f.contig.compartment).tab().append(view)
			.tab().append(role).tab().append(f.id).tab().append(f.label).tab().append(f.start).tab().append(f.stop).tab().append(f.strand).tab().append(outcome).tab().append(partner==null?"NA":partner.id);
		if(partner==null || !outcome.equals("TP")){b.append("\tNA\tNA");}else{b.tab().append(d5).tab().append(d3);}
		final Genome g=f.contig.genome;final String kind=g.srpReferenceKind!=null && view.startsWith("SRP")?g.srpReferenceKind:g.referenceKind;
		out.print(b.tab().append(f.geometry).tab().append(g.panel.equals("P1")?"SOURCE_ANNOTATION_RELATIVE":kind).nl());
	}
	static void unknown(ByteStreamWriter out,ByteBuilder b,String arm,Feature f,String role,String reason){
		assert(f!=null):"Unknown classification is retained with original coordinates and attributes";
		out.print(b.clear().append(arm).tab().append(f.contig.genome.panel).tab().append(f.contig.genome.id).tab().append(f.contig.id).tab().append(role).tab().append(f.id)
			.tab().append(f.start).tab().append(f.stop).tab().append(f.strand).tab().append(reason).tab().append(f.attributes).nl());
	}
	static void readTime(String path,ByteStreamWriter out,ByteBuilder b,String arm,String panel,String task,int genomes){
		assert(genomes>0):"Runtime is per process; parallel array wall times must not be presented as panel elapsed time";
		final ByteFile in=RrnaResourceIO.open(path);String line;try{final byte[] raw=in.nextLine();require(raw!=null && in.nextLine()==null,"One process timing row required");line=new String(raw,java.nio.charset.StandardCharsets.US_ASCII);}finally{require(!in.close(),"Timing read failed");}
		final String[] tokens=line.split(" ");require(tokens.length==3 && tokens[0].startsWith("wall_seconds=") && tokens[1].startsWith("user_seconds=") && tokens[2].startsWith("system_seconds="),"Expected literal time producer schema");
		b.clear().append(arm).tab().append(panel).tab().append(task).tab().append(genomes);for(String t:tokens){final double value=Double.parseDouble(t.substring(t.indexOf('=')+1));require(Double.isFinite(value) && value>=0,"Invalid process time");b.tab().append(value,6);}out.print(b.nl());
	}
	static int delta(int d,int length){if(length<=0){return d;}int x=d%length;if(x>length/2){x-=length;}if(x< -length/2){x+=length;}return x;}
	public static void readWorkload(String path,ByteStreamWriter out,ByteBuilder b,String arm,String panel,String task){
		readWorkload(path,out,b,arm,panel,task,true);
	}
	/** Ordinary CallGenes logs omit tRNA alignments unless extendedStats/verbose
	 * is enabled (CallGenes.printStats); optional means unavailable, never zero. */
	public static void readWorkload(String path,ByteStreamWriter out,ByteBuilder b,String arm,String panel,String task,boolean requireTrna){
		assert(out!=null):"Workload belongs beside timings for each saved caller process, not a new caller run";
		final ByteFile in=RrnaResourceIO.open(path);final ObjectSet<String> seen=new ObjectSet<String>(String.class);int families=0;
		try{for(byte[] line;(line=in.nextLine())!=null;){
			String metric=null,family=null;int colon=-1;
			if(Tools.startsWith(line,"ncRNA Alignments (")){
				final int start="ncRNA Alignments (".length();int end=start;while(end<line.length && line[end]!=')'){end++;}
				require(end>start && end+1<line.length && line[end+1]==':',"Malformed family alignment counter in "+path);
				family=new String(line,start,end-start,java.nio.charset.StandardCharsets.US_ASCII);metric="alignments";colon=end+1;
				require(known(fold(family)),"Unknown workload family: "+family);families++;
			}else if(Tools.startsWith(line,"tRNA Alignments:")){family="tRNA";metric="alignments";colon="tRNA Alignments".length();}
			else if(Tools.startsWith(line,"ncRNA Shared Sweep Passes:")){family="GLOBAL";metric="shared_sweep_passes";colon="ncRNA Shared Sweep Passes".length();}
			else if(Tools.startsWith(line,"ncRNA Aligned Cells:")){family="GLOBAL";metric="aligned_cells";colon="ncRNA Aligned Cells".length();}
			else{continue;}
			// The colon offset is derived from the literal prefix, then whitespace is trimmed.
			require(colon>=0 && line[colon]==':',"Missing workload counter separator");
			int from=colon+1,to=line.length;while(from<to && line[from]<=32){from++;}while(to>from && line[to-1]<=32){to--;}
			require(from<to && to-from<=18,"Missing or oversized workload counter value");for(int i=from;i<to;i++){require(line[i]>='0' && line[i]<='9',"Workload counter must be an unsigned integer");}
			final long count=parse.Parse.parseLong(line,from,to);require(count>=0 && seen.add(metric+'\t'+family),"Negative/duplicate workload counter");
			out.print(b.clear().append(arm).tab().append(panel).tab().append(task).tab().append(metric).tab().append(family).tab().append(count).nl());
		}}finally{require(!in.close(),"Caller workload log read failed");}
		final boolean hasTrna=seen.contains("alignments\ttRNA");
		require(families>0 && (!requireTrna || hasTrna) && seen.contains("shared_sweep_passes\tGLOBAL") && seen.contains("aligned_cells\tGLOBAL"),"Incomplete caller workload counters in "+path);
		if(!hasTrna){out.print(b.clear().append(arm).tab().append(panel).tab().append(task).append("\talignments\ttRNA\tNA\n"));}
	}
	static void fraction(ByteBuilder b,long n,long d){b.tab();if(d==0){b.append("NA");}else{b.append(n/(double)d,8);}}
	static boolean valid(Contig c,int start,int stop,boolean source){return start>0 && start<=c.length && stop>=start && (long)stop-start+1<=c.length && (source || stop<=c.length);}
	static boolean rna(LineParser1 p){
		assert(p.terms()==9):"RNA filtering reads only the original GFF type field before allocating retained identities";
		if(p.termEquals("mRNA",2)){return false;}p.setBounds(2);final byte[] line=p.line();final int end=p.b();
		return end-p.a()>=3 && line[end-3]=='R' && line[end-2]=='N' && line[end-1]=='A';
	}
	static String fold(String label){
		if(label.equals("16S") || label.equals("18S") || label.equals("euk18S")){return "SSU";}
		if(label.equals("23S") || label.equals("lsu") || label.equals("eukLSU")){return "LSU";}
		if(label.equals("prok5S") || label.equals("euk5S")){return "5S";}if(label.equals("r58")){return "5.8S";}
		if(label.equals("rnasep")){return "RNaseP";}if(label.equals("tmrna")){return "tmRNA";}
		if(label.equals("srp_small")){return "SRP_small";}if(label.equals("srp_large")){return "SRP_large";}
		if(label.equals("sixs_rf00013") || label.equals("sixs_rf01685")){return "6S";}
		return known(label)?label:"UNKNOWN";
	}
	static boolean known(String family){if(sixsView(family)){return true;}for(String f:FAMILIES){if(f.equals(family)){return true;}}return false;}
	static boolean sixsView(String view){return view.equals("6S_RF00013") || view.equals("6S_RF01685");}
	public static boolean inView(Feature f,String view){
		assert(f!=null):"Subtype views use the retained CM accession/caller label, never length guesses";
		if(sixsView(view)){return f.family.equals("6S") && (f.label.equals(view.substring(3)) || f.label.equals("sixs_"+view.substring(3).toLowerCase(Locale.ROOT)));}
		return inView(f.family,view);
	}
	static boolean inView(String family,String view){return family.equals(view) || view.equals("SRP") && (family.equals("SRP_small") || family.equals("SRP_large"));}
	static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	public static final String[] FAMILIES={"5S","5.8S","SSU","LSU","tRNA","RNaseP","tmRNA","6S","SRP_small","SRP_large","SRP"};
	public static final String HISTORY_HEADER="assembly\tnative_seqid\tstrand\tstart\tend\tfold_family\thistory_class\tmask_eligible\tsource_gff\tsource_sha80\tsource_line\toriginal_type\toriginal_attributes\tmask_reason";
	public static final class Genome{
		public Genome(String i,String p,String f,String h){id=i;panel=p;fasta=f;hash=h;}
		final String id,panel,fasta,hash;public String source,sourceHash;
		String referenceKind="RAW_GA_HITS_PROVISIONAL";
		String srpReferenceKind;
		public final ObjectMap<String,Contig> contigs=new ObjectMap<String,Contig>(String.class,Contig.class);public final ArrayList<Contig> ordered=new ArrayList<Contig>();
	}
	public static final class Contig{
		public Contig(Genome g,String i,int n,String c){require(n>0,"Positive contig length required");genome=g;id=i;length=n;compartment=c;}
		public final Genome genome;public final String id,compartment;final int length;public boolean circular;
		public final ArrayList<Feature> truth=new ArrayList<Feature>(),masks=new ArrayList<Feature>(),unknown=new ArrayList<Feature>(),calls=new ArrayList<Feature>();
	}
	public static final class Feature{
		public Feature(Contig c,String f,String label_,int a,int z,String s,String i,String at,String origin_){contig=c;family=f;label=label_;start=a;stop=z;strand=s;id=i;attributes=at;origin=origin_;}
		public Locus locus(String view){return new Locus(contig.id,view,start,stop,strand.charAt(0));}
		final Contig contig;public final String family,label,strand,id,attributes,origin;final int start,stop;public String geometry="SOURCE_UNSTRATIFIED",state="CANDIDATE";public boolean truncated;
		boolean boundaryEligible=true,bothEligible=true;
	}
	public final ObjectMap<String,String> modelFamilies=new ObjectMap<String,String>(String.class,String.class);
	public final ObjectMap<String,Genome> genomes=new ObjectMap<String,Genome>(String.class,Genome.class);final ArrayList<Genome> ordered=new ArrayList<Genome>();long rawHits;
	public final ObjectMap<String,long[]> totals=new ObjectMap<String,long[]>(String.class,long[].class),truthStrata=new ObjectMap<String,long[]>(String.class,long[].class),labelCounts=new ObjectMap<String,long[]>(String.class,long[].class);
	public final ObjectMap<String,long[]> physicalCounts=new ObjectMap<String,long[]>(String.class,long[].class);
	final ObjectMap<String,long[]> labelUniverse=new ObjectMap<String,long[]>(String.class,long[].class),armPanels=new ObjectMap<String,long[]>(String.class,long[].class);
	final ObjectSet<String> covered=new ObjectSet<String>(String.class);final ObjectMap<String,long[]> coverageCounts=new ObjectMap<String,long[]>(String.class,long[].class);
	String membershipPath;
	boolean frozenMode,historyLoaded;long expectedFrozenLoci,referenceRows,historyRows,historicalMasks;
	String frozenReferenceKind="FROZEN_LOCUS_TRUTH_V1";
	final ObjectSet<String> membershipMembers=new ObjectSet<String>(String.class),membershipArms=new ObjectSet<String>(String.class);
	final ObjectMap<String,long[]> membershipPanels=new ObjectMap<String,long[]>(String.class,long[].class);
	static final Comparator<Feature> CALL_ORDER=(a,b)->{int c=a.strand.compareTo(b.strand);if(c==0){c=Integer.compare(a.start,b.start);}if(c==0){c=Integer.compare(a.stop,b.stop);}if(c==0){c=a.id.compareTo(b.id);}return c;};
}
