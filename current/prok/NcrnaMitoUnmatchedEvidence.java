package prok;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Locale;
import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import stream.Writer;
import stream.WriterFactory;
import structures.ByteBuilder;
import tax.TaxTree;

/** Saved-call evidence, without new calling, alignment, truth editing or origin inference.
 * Selects newly unmatched nuclear SSU calls and retains the complete mitochondrial
 * truth cohort. Every overlap keeps its original row and strand relation.
 * The committed producer verifies source hashes; this helper retains their
 * original roster bindings and does not substitute a new biological truth.
 * @author Brian, Ganyu */
public final class NcrnaMitoUnmatchedEvidence {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args,NcrnaMitoUnmatchedEvidence.class,false);
		final ObjectMap<String,String> flags=new ObjectMap<String,String>(String.class,String.class);
		try{
			for(String arg:pp.args){
				final String[] parts=Parse.splitOnFirst(arg,'=');
				require(parts.length==2,"Expected flag=value");final String key=parts[0].toLowerCase(Locale.ROOT);
				require(KNOWN.contains("|"+key+"|") && flags.put(key,parts[1])==null,"Unknown/repeated argument: "+key);
			}
			for(String key:REQUIRED){require(flags.get(key)!=null,"Missing evidence binding: "+key);}
			final int expected=Integer.parseInt(flags.get("expected")),mito=Integer.parseInt(flags.get("expectedmito"));
			require(expected>0 && mito>0,"Both complete reviewed populations must have positive denominators");
			Shared.setThreads(1);ReadWrite.setZipThreads(1);Shared.TRIM_READ_DESCRIPTION=false;Read.TO_UPPER_CASE=false;
			new NcrnaMitoUnmatchedEvidence().run(flags,expected,mito);
		}finally{Shared.closeStream(pp.outstream);}
	}

	public void run(ObjectMap<String,String> flags,int expected,int expectedMito){
		assert(cases.isEmpty()):"One immutable comparison per helper instance";
		Shared.setThreads(1);ReadWrite.setZipThreads(1);Shared.TRIM_READ_DESCRIPTION=false;Read.TO_UPPER_CASE=false;
		final String prefix=flags.get("outprefix");
		RrnaResourceIO.fresh(prefix,".cases.tsv",".overlaps.tsv.gz",".contigs.tsv",".contigs.fa.gz",".spans.fa.gz",".summary.tsv");
		loadAudit(flags.get("audit"),expected,expectedMito);
		loadInputs(flags.get("inputs"));
		bindVerifiedInputs(flags.get("verified"));
		loadTaxonomy(flags.get("tree"));
		loadModels(flags.get("models"));
		final ByteStreamWriter overlaps=RrnaResourceIO.writer(prefix+".overlaps.tsv.gz");
		try{
			overlaps.println("case_id\tevidence_kind\tsource_path\tsource_line\toverlap_bp\tsame_strand\toriginal_row");
			for(String panel:new String[]{"p2","p3"}){for(int arm=0;arm<ARMS.length;arm++){
				loadCalls(flags.get(panel+ARMS[arm].toLowerCase(Locale.ROOT)),panel.equals("p2")?"P2":"P3",arm,overlaps);
			}}
			loadTruth(flags.get("truth"),overlaps);
			loadHistory(flags.get("history"),overlaps);
			loadMasks(flags.get("masks"),overlaps);
		}finally{RrnaResourceIO.finish(overlaps);}
		loadSequences(prefix);
		write(prefix,expected,expectedMito);
	}

	private void loadAudit(String path,int expected,int expectedMito){
		final ObjectSet<String> old=new ObjectSet<String>(String.class);
		final ArrayList<Case> candidates=new ArrayList<Case>();
		final ObjectMap<String,String[]> mitoStates=new ObjectMap<String,String[]>(String.class,String[].class);
		final ObjectMap<String,Case> mitoRows=new ObjectMap<String,Case>(String.class,Case.class);
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{
			RrnaResourceIO.header(in,NcrnaP1MissInventory.AUDIT_HEADER,path);
			for(byte[] row;(row=in.nextLine())!=null;){
				p.set(row);require(p.terms()==18,"Native saved audit must retain all18 columns");
				if(!p.termEquals("SSU",5)){continue;}
				if(p.termEquals("P2",1) && p.termEquals("NUCLEAR",4) && p.termEquals("CALL",6)){
					if(p.termEquals("baseline",0)){old.add(coordinate(p));}
					else if(p.termEquals("B",0) && (p.termEquals("APPARENT_FP",12) || p.termEquals("DUPLICATE_APPARENT_FP",12))){candidates.add(new Case(p,"NUCLEAR_NEW_UNMATCHED"));}
				}else if(p.termEquals("P3",1) && p.termEquals("MITOCHONDRION",4) && p.termEquals("TRUTH",6)){
					final String id=p.parseString(7);final int arm=arm(p.parseString(0));
					String[] states=mitoStates.get(id);if(states==null){states=new String[3];mitoStates.put(id,states);}
					require(states[arm]==null && (p.termEquals("TP",12) || p.termEquals("FN",12)),"Every mitochondrial truth has exactly one outcome per arm");
					states[arm]=p.parseString(12);
					final Case c=new Case(p,"ORGANELLE_SSU_TRUTH"),prior=mitoRows.get(id);
					if(prior==null){mitoRows.put(id,c);}else{require(prior.coordinate.equals(c.coordinate),"Mitochondrial reference coordinates cannot change across arms");}
				}
			}
		}finally{require(!in.close(),"Saved audit read failed");}
		for(Case c:candidates){if(!old.contains(c.coordinate)){cases.add(c);}}
		require(cases.size()==expected,"New unmatched coordinate population changed: "+cases.size()+" versus "+expected);
		for(String id:sorted(mitoRows)){
			final Case c=mitoRows.get(id);c.states=mitoStates.get(id);
			for(String state:c.states){require(state!=null,"All three arms must report every mitochondrial truth, including FN");}cases.add(c);
		}
		require(mitoRows.size()==expectedMito,"Complete mitochondrial truth population changed");
		cases.sort((a,b)->{int n=a.kind.compareTo(b.kind);if(n==0){n=a.assembly.compareTo(b.assembly);}if(n==0){n=a.contig.compareTo(b.contig);}if(n==0){n=Integer.compare(a.start,b.start);}if(n==0){n=Integer.compare(a.stop,b.stop);}if(n==0){n=a.strand.compareTo(b.strand);}return n;});
		final ObjectSet<String> unique=new ObjectSet<String>(String.class);
		for(int i=0;i<cases.size();i++){
			final Case c=cases.get(i);c.id="MC"+(i+1);require(unique.add(c.kind+'\t'+c.coordinate),"Duplicate target coordinates require explicit multiplicity handling");
			final String key=c.assembly+'\t'+c.contig;TargetContig t=contigs.get(key);
			if(t==null){t=new TargetContig(c.assembly,c.contig,"MCONTIG_"+(contigs.size()+1));contigs.put(key,t);}t.cases.add(c);c.source=t;
			if(c.kind.equals("NUCLEAR_NEW_UNMATCHED")){
				t.export=true;final int colon=c.auditId.lastIndexOf(':');require(colon>0,"Native call ID must preserve source path and physical line");
				c.callPath=c.auditId.substring(0,colon);c.callLine=Long.parseLong(c.auditId.substring(colon+1));
				require(c.callLine>0 && callLines.get(c.callLine)==null,"One candidate GFF line per selected call");callLines.put(c.callLine,c);
			}
		}
	}

	private void loadInputs(String path){
		assert(!contigs.isEmpty()):"Read only assemblies implicated by native evidence";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{
			p.set(in.nextLine());require(p.terms()==15 && p.termEquals("assembly",1) && p.termEquals("taxid",3) && p.termEquals("fasta",6),"Expected original198-panel input schema; host taxid is the fourth column");
			for(byte[] row;(row=in.nextLine())!=null;){p.set(row);require(p.terms()==15,"Panel input schema changed");
				final String id=p.parseString(1);ArrayList<TargetContig> selected=null;
				for(TargetContig t:contigs.values()){if(t!=null && t.assembly.equals(id)){if(selected==null){selected=new ArrayList<TargetContig>();}selected.add(t);}}
				if(selected!=null){final Assembly a=new Assembly(id,p.parseString(6),p.parseString(7),p.parseInt(3),selected,path);require(assemblies.put(id,a)==null,"Selected source assembly must be unique: "+id);for(TargetContig t:selected){t.input=a;}}
			}
		}finally{require(!in.close(),"Original input roster read failed");}
		for(TargetContig t:contigs.values()){if(t!=null){require(t.input!=null,"Missing source FASTA for "+t.assembly);}}
	}

	private void bindVerifiedInputs(String path){
		require(path!=null,"Bind the accepted original-source hash ledger before resolving roster NA hashes");
		final ObjectMap<String,String> hashes=EukPanelGffValidation.hashes(path);
		for(Assembly a:assemblies.values()){if(a!=null){
			final String hash=hashes.get(a.fasta);
			require(hash!=null,"No verified hash for exact FASTA path: assembly="+a.id+" fasta="+a.fasta+" taxid="+a.taxid+" rosterHash="+a.originalHash);
			require(a.originalHash.equals("NA") || a.originalHash.equals(hash),"Verified FASTA hash disagrees with original roster: assembly="+a.id+" roster="+a.originalHash+" verified="+hash);
			a.hash=hash;a.hashSource=path;
		}}
	}

	private void loadTaxonomy(String path){
		require(path!=null && !path.equalsIgnoreCase("auto"),"Host lineage requires the producer's pinned taxonomy tree, not an implicit default");
		final TaxTree tree;final boolean old=TaxTree.LOAD_FAKE_TREE;
		try{TaxTree.LOAD_FAKE_TREE=false;tree=TaxTree.loadTaxTreeFile(path);}finally{TaxTree.LOAD_FAKE_TREE=old;}
		for(Assembly a:assemblies.values()){if(a!=null){
			a.phylum=NcrnaMitoClusterEvaluation.lineage(tree,a.taxid,TaxTree.PHYLUM_E);
			a.hostClass=NcrnaMitoClusterEvaluation.lineage(tree,a.taxid,TaxTree.CLASS_E);
		}}
	}

	private void loadModels(String path){
		RrnaResourceIO.stream(path,r->{final String name=RrnaResourceIO.first(r.id);require(modelLengths.put(name,new int[]{r.length()})==null,"Winning model names must be unique");});
		require(modelLengths.size()>0,"Bind the actual candidate16S consensus library");
	}

	private void loadCalls(String path,String panel,int arm,ByteStreamWriter out){
		assert(path!=null && arm>=0 && arm<3):"Overlap rows refer to original caller outputs for each arm";
		final ObjectMap<String,ArrayList<Case>> selected=new ObjectMap(String.class,ArrayList.class);
		for(Case c:cases){if(c.panel.equals(panel)){
			if(arm==2 && c.kind.equals("NUCLEAR_NEW_UNMATCHED")){require(path.equals(c.callPath),"Selected call IDs must bind this exact original candidate GFF");}
			ArrayList<Case> list=selected.get(c.contig);if(list==null){list=new ArrayList<Case>();selected.put(c.contig,list);}
			else{require(list.get(0).assembly.equals(c.assembly),"A batched GFF cannot resolve repeated contig IDs across assemblies");}list.add(c);
		}}
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');long line=0;
		String currentSeq=null;ArrayList<Case> currentTargets=null;
		try{for(byte[] row;(row=in.nextLine())!=null;){line++;if(row.length==0 || row[0]=='#'){continue;}p.set(row);require(p.terms()==9,"Saved calls must retain fullGFF9 rows");
			final Case exact=panel.equals("P2") && arm==2 ? callLines.get(line) : null;
			if(exact!=null){
				require(panel.equals("P2") && arm==2 && p.termEquals(exact.contig,0) && p.parseInt(3)==exact.start && p.parseInt(4)==exact.stop && p.termEquals(exact.strand,6),"Audit call ID must bind its exact unmodified GFF row");
				require(exact.model==null,"One physical GFF row per selected call");exact.attributes=p.parseString(8);exact.model=RrnaShortlistReplay.model(exact.attributes);
				require("16S".equals(PanelAnnotationGrade.colonField(exact.attributes,"family")) && modelLengths.get(exact.model)!=null,"Selected added call must name a loaded16S model");
				exact.score=p.parseString(5);exact.rawGff=new String(row,java.nio.charset.StandardCharsets.UTF_8);
			}
			if(currentSeq==null || !p.termEquals(currentSeq,0)){currentSeq=p.parseString(0);currentTargets=selected.get(currentSeq);}
			if(currentTargets==null){continue;}final boolean rna=NcrnaClippingComparison.rna(p),cds=p.termEquals("CDS",2);
			for(Case c:currentTargets){final int n=overlap(c.start,c.stop,p.parseInt(3),p.parseInt(4));if(n>0){
				if(rna){c.rnaOverlaps[arm]++;}else if(cds){c.cdsOverlaps[arm]++;}else{c.otherOverlaps[arm]++;}
				evidence(out,c,(rna?"CALL_RNA_":cds?"CALL_CDS_":"CALL_OTHER_")+ARMS[arm],path,line,n,p.termEquals(c.strand,6),row);
			}}
		}}finally{require(!in.close(),"Saved full caller GFF read failed");}
	}

	private void loadTruth(String path,ByteStreamWriter out){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');long line=1;
		try{RrnaResourceIO.header(in,EukPanelTruthFreeze.HEADER,path);for(byte[] row;(row=in.nextLine())!=null;){line++;p.set(row);
			require(p.terms()==35 && p.termEquals("TRUTH_V2",32),"Only frozen truth_v2 supplies current reference evidence");
			final TargetContig t=contigs.get(p.parseString(1)+'\t'+p.parseString(2));if(t==null){continue;}
			for(Case c:t.cases){final int n=overlap(c.start,c.stop,p.parseInt(4),p.parseInt(5));if(n>0){
				final boolean admitted=p.termEquals("ADMITTED",26);require(admitted || p.termEquals("BELOW_GA",26),"Frozen admission must not silently reinterpret an unknown label");if(admitted){c.truthOverlaps++;}else{c.maskOverlaps++;}
				evidence(out,c,admitted?"TRUTH_ADMITTED":"TRUTH_BELOW_GA",path,line,n,p.termEquals(c.strand,3),row);
			}}
		}}finally{require(!in.close(),"Frozen truth read failed");}
	}
	private void loadHistory(String path,ByteStreamWriter out){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');long line=1;
		try{RrnaResourceIO.header(in,NcrnaClippingComparison.HISTORY_HEADER,path);for(byte[] row;(row=in.nextLine())!=null;){line++;p.set(row);require(p.terms()==14,"Historical mask schema changed");
			final TargetContig t=contigs.get(p.parseString(0)+'\t'+p.parseString(1));if(t==null){continue;}
			for(Case c:t.cases){final int n=overlap(c.start,c.stop,p.parseInt(3),p.parseInt(4));if(n>0){
				final boolean eligible=EukPanelLocusDiagnostics.bit(p,7);if(eligible){c.maskOverlaps++;}c.historyOverlaps++;
				evidence(out,c,eligible?"HISTORICAL_MASK":"HISTORICAL_NONMASK",path,line,n,p.termEquals(c.strand,2),row);
			}}
		}}finally{require(!in.close(),"Historical evidence read failed");}
	}
	private void loadMasks(String path,ByteStreamWriter out){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');long line=2;
		try{final byte[] comment=in.nextLine();require(comment!=null && comment.length>0 && comment[0]=='#',"Expected source-mask policy comment");p.set(in.nextLine());require(p.terms()==33 && p.termEquals("mask_eligible",31),"Expected source-mask schema");
			for(byte[] row;(row=in.nextLine())!=null;){line++;p.set(row);require(p.terms()==33,"Source-mask schema changed");
				final TargetContig t=contigs.get(p.parseString(1)+'\t'+p.parseString(8));if(t==null){continue;}
				for(Case c:t.cases){final int n=overlap(c.start,c.stop,p.parseInt(17),p.parseInt(18));if(n>0){
					final boolean eligible=EukPanelLocusDiagnostics.bit(p,31);if(eligible){c.maskOverlaps++;}
					evidence(out,c,eligible?"OUT_OF_MAP_MASK":"OUT_OF_MAP_NONMASK",path,line,n,p.termEquals(c.strand,20),row);
				}}
			}
		}finally{require(!in.close(),"Outside-map evidence read failed");}
	}

	private void loadSequences(String prefix){
		final Writer whole=WriterFactory.makeWriter(FileFormat.testOutput(prefix+".contigs.fa.gz",FileFormat.FASTA,null,true,false,false,true),null,1);
		final Writer spans=WriterFactory.makeWriter(FileFormat.testOutput(prefix+".spans.fa.gz",FileFormat.FASTA,null,true,false,false,true),null,1);
		whole.start();spans.start();final long[] batches={0,0};
		try{for(String id:sorted(assemblies)){
			final Assembly a=assemblies.get(id);final ObjectMap<String,TargetContig> wanted=new ObjectMap<String,TargetContig>(String.class,TargetContig.class);
			for(TargetContig t:a.contigs){require(wanted.put(t.seqid,t)==null,"Target contigs must be unique within an assembly");}
			RrnaResourceIO.stream(a.fasta,r->{final TargetContig t=wanted.get(RrnaResourceIO.first(r.id));if(t==null){return;}
				require(t.header==null && r.bases!=null,"Each selected source contig must occur exactly once");t.header=r.id;t.length=r.length();
				for(byte base:r.bases){final int n=RrnaResourceIO.base(base);if(n>=0){t.acgt++;if(n==1 || n==2){t.gc++;}}}
				if(t.export){final ArrayList<Read> batch=new ArrayList<Read>(1);batch.add(new Read(r.bases,null,t.exportId+" assembly="+t.assembly+" original_header="+r.id,batches[0],false));whole.add(batch,batches[0]++);}
				final ArrayList<Read> batch=new ArrayList<Read>(t.cases.size());
				for(Case c:t.cases){require(c.start>=1 && c.stop>=c.start && c.stop<=r.length(),"Call/truth coordinates must fit the actual original FASTA");
					final byte[] seq=java.util.Arrays.copyOfRange(r.bases,c.start-1,c.stop);if(c.strand.equals("-")){AminoAcid.reverseComplementBasesInPlace(seq);}
					batch.add(new Read(seq,null,c.id+" assembly="+c.assembly+" contig="+c.contig+" start="+c.start+" stop="+c.stop+" strand="+c.strand,batches[1],false));}
				spans.add(batch,batches[1]++);
			});
		}}finally{final boolean failed=whole.poisonAndWait()|spans.poisonAndWait();require(!failed,"Complete native FASTA exports are required before evidence is reported");}
		for(TargetContig t:contigs.values()){if(t!=null){require(t.header!=null,"Selected contig missing from original FASTA: "+t.assembly+"/"+t.seqid);}}
	}

	private void write(String prefix,int expected,int mito){
		final ByteStreamWriter out=RrnaResourceIO.writer(prefix+".cases.tsv"),ct=RrnaResourceIO.writer(prefix+".contigs.tsv"),summary=RrnaResourceIO.writer(prefix+".summary.tsv");
		final ByteBuilder b=new ByteBuilder();int exported=0;long bases=0;
		try{
			out.println("case_id\tkind\tassembly\tcontig\tstart1\tstop1\tstrand\tlength\tmodel\tGFF_score\tverifier_identity\tidentity_provenance\tcontig_length\tcontig_ACGT\tcontig_GC\tGC_fraction_ACGT\ttruth_overlaps_any_strand\tmask_overlaps_any_strand\thistorical_overlaps_any_strand\tbaseline_RNA_overlaps\tA_RNA_overlaps\tB_RNA_overlaps\tbaseline_truth_outcome\tA_truth_outcome\tB_truth_outcome\torigin\tsource_fasta\tsource_sha80\toriginal_header\taudit_id\toriginal_GFF9\tbaseline_CDS_overlaps\tA_CDS_overlaps\tB_CDS_overlaps\tbaseline_other_overlaps\tA_other_overlaps\tB_other_overlaps\thost_taxid\thost_phylum\thost_class\toriginal_roster_sha80\tsource_hash_binding");
			for(Case c:cases){final TargetContig t=c.source;
				if(c.kind.equals("NUCLEAR_NEW_UNMATCHED")){require(c.model!=null,"Every selected audit call must bind its original GFF record");}
				b.clear().append(c.id).tab().append(c.kind).tab().append(c.assembly).tab().append(c.contig).tab().append(c.start).tab().append(c.stop).tab().append(c.strand).tab().append(c.stop-c.start+1)
					.tab().append(value(c.model)).tab().append(value(c.score)).append("\tNA\tORIGINAL_VERIFIER_NOT_RECORDED\t").append(t.length).tab().append(t.acgt).tab().append(t.gc).tab();
				fraction(b,t.gc,t.acgt);b.tab().append(c.truthOverlaps).tab().append(c.maskOverlaps).tab().append(c.historyOverlaps);
				for(int n:c.rnaOverlaps){b.tab().append(n);}for(int i=0;i<3;i++){b.tab().append(c.states==null?"NA":c.states[i]);}
				b.append("\tUNRESOLVED\t").append(t.input.fasta).tab().append(t.input.hash).tab().append(escape(t.header)).tab().append(c.auditId).tab().append(escape(value(c.rawGff)));
				for(int n:c.cdsOverlaps){b.tab().append(n);}for(int n:c.otherOverlaps){b.tab().append(n);}
				out.print(b.tab().append(t.input.taxid).tab().append(t.input.phylum).tab().append(t.input.hostClass).tab().append(t.input.originalHash).tab().append(escape(t.input.hashSource)).nl());
			}
			ct.println("assembly\tcontig\tlength\tACGT\tGC\tGC_fraction_ACGT\texported_id\tsource_fasta\tsource_sha80\toriginal_header\toriginal_roster_sha80\tsource_hash_binding");
			for(String key:sorted(contigs)){final TargetContig t=contigs.get(key);b.clear().append(t.assembly).tab().append(t.seqid).tab().append(t.length).tab().append(t.acgt).tab().append(t.gc).tab();fraction(b,t.gc,t.acgt);
				ct.print(b.tab().append(t.export?t.exportId:"NA").tab().append(t.input.fasta).tab().append(t.input.hash).tab().append(escape(t.header)).tab().append(t.input.originalHash).tab().append(escape(t.input.hashSource)).nl());if(t.export){exported++;bases+=t.length;}}
			summary.println("new_unmatched_calls\torganelle_truth_loci\tsource_assemblies\texported_full_contigs\texported_bases\tnew_alignments\ttruth_edits\torigin_classification");
			summary.print(b.clear().append(expected).tab().append(mito).tab().append(assemblies.size()).tab().append(exported).tab().append(bases).append("\t0\t0\tUNRESOLVED\n"));
		}finally{RrnaResourceIO.finish(out,ct,summary);}
		System.err.println("MITO_UNMATCHED_EVIDENCE_PASS added="+expected+" organelleTruth="+mito+" fullContigs="+exported+" no_alignment_no_truth_edit=true");
	}
	private void evidence(ByteStreamWriter out,Case c,String kind,String path,long line,int overlap,boolean same,byte[] row){
		assert(overlap>0):"Only actual geometric overlap belongs in the retained evidence rows";
		out.print(evidenceBuffer.clear().append(c.id).tab().append(kind).tab().append(path).tab().append(line).tab().append(overlap).tab().append(same).tab().append(row).nl());
	}
	private static <V> ArrayList<String> sorted(ObjectMap<String,V> map){
		assert(map!=null):"Deterministic output order cannot depend on hash iteration";
		final ArrayList<String> keys=new ArrayList<String>();for(String key:map.keys()){if(key!=null){keys.add(key);}}Collections.sort(keys);return keys;
	}
	public static int overlap(int a,int z,int b,int y){require(a>=1 && z>=a && b>=1 && y>=b,"Valid inclusive source intervals required");return Math.max(0,Math.min(z,y)-Math.max(a,b)+1);}
	private static String coordinate(LineParser1 p){return p.parseString(2)+'\t'+p.parseString(3)+'\t'+p.parseInt(9)+'\t'+p.parseInt(10)+'\t'+p.parseString(11);}
	private static int arm(String s){for(int i=0;i<ARMS.length;i++){if(ARMS[i].equals(s)){return i;}}throw new IllegalArgumentException("Unexpected saved arm: "+s);}
	private static String value(String s){return s==null?"NA":s;}
	private static String escape(String s){return s.replace("\\","\\\\").replace("\t","\\t").replace("\r","\\r").replace("\n","\\n");}
	private static void fraction(ByteBuilder b,long n,long d){if(d==0){b.append("NA");}else{b.append(n/(double)d,9);}}
	private static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final String[] ARMS={"baseline","A","B"};
	private static final String[] REQUIRED={"audit","inputs","verified","models","tree","truth","history","masks","p2baseline","p2a","p2b","p3baseline","p3a","p3b","outprefix","expected","expectedmito"};
	private static final String KNOWN="|audit|inputs|verified|models|tree|truth|history|masks|p2baseline|p2a|p2b|p3baseline|p3a|p3b|outprefix|expected|expectedmito|";
	private static final class Case {
		Case(LineParser1 p,String kind_){kind=kind_;panel=p.parseString(1);assembly=p.parseString(2);contig=p.parseString(3);auditId=p.parseString(7);start=p.parseInt(9);stop=p.parseInt(10);strand=p.parseString(11);coordinate=coordinate(p);require(start>0 && stop>=start && (strand.equals("+") || strand.equals("-")),"Native case geometry must be valid");}
		final String kind,panel,assembly,contig,auditId,strand,coordinate;final int start,stop;String id,model,score,attributes,rawGff,callPath;long callLine;String[] states;TargetContig source;int truthOverlaps,maskOverlaps,historyOverlaps;final int[] rnaOverlaps=new int[3],cdsOverlaps=new int[3],otherOverlaps=new int[3];
	}
	private static final class TargetContig {
		TargetContig(String a,String s,String e){assembly=a;seqid=s;exportId=e;}
		final String assembly,seqid,exportId;final ArrayList<Case> cases=new ArrayList<Case>();Assembly input;String header;int length;long acgt,gc;boolean export;
	}
	private static final class Assembly {
		Assembly(String i,String f,String h,int tid,ArrayList<TargetContig> c,String roster){
			id=i;fasta=f;hash=originalHash=h;taxid=tid;contigs=c;hashSource=roster;organelle=c.get(0).cases.get(0).panel.equals("P3");
			require(tid>0,"Missing positive host taxid in source roster: assembly="+id+" taxid="+tid+" hash="+h);
			require(SeedOffsetTable.validSha80(h) || (organelle && "NA".equals(h)),"Invalid original FASTA hash: assembly="+id+" taxid="+tid+" hash="+h+"; only P3 roster NA may resolve against the verified source ledger");
		}
		final String id,fasta,originalHash;final int taxid;final boolean organelle;final ArrayList<TargetContig> contigs;String hash,hashSource,phylum,hostClass;
	}
	private final ArrayList<Case> cases=new ArrayList<Case>();
	private final map.LongObjectMap<Case> callLines=new map.LongObjectMap<Case>(128,Case.class);
	private final ByteBuilder evidenceBuffer=new ByteBuilder();
	private final ObjectMap<String,TargetContig> contigs=new ObjectMap<String,TargetContig>(String.class,TargetContig.class);
	private final ObjectMap<String,Assembly> assemblies=new ObjectMap<String,Assembly>(String.class,Assembly.class);
	private final ObjectMap<String,int[]> modelLengths=new ObjectMap<String,int[]>(String.class,int[].class);
}
