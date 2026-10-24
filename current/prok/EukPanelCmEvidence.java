package prok;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Retains independent CM evidence and physical historical overlaps, never admits truth.
 * The producer verifies hashes with canonical sha80 before supplying verified=.
 * Model identity includes the CM bytes; equal accessions alone do not establish equality.
 * @author Ganyu, Raiden */
public final class EukPanelCmEvidence {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args, EukPanelCmEvidence.class, false);
		String inputs=null, sources=null, contigs=null, contigSources=null, features=null;
		String runs=null, models=null, verified=null, out=null;int expected=198, expectedModels=25;boolean prokInputs=false;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0,"Expected flag=value");
			final String k=arg.substring(0,eq).toLowerCase(Locale.ROOT), v=arg.substring(eq+1);
			if(k.equals("inputs")){inputs=v;}else if(k.equals("sources")){sources=v;}
			else if(k.equals("contigs")){contigs=v;}else if(k.equals("contigsources")){contigSources=v;}
			else if(k.equals("features")){features=v;}else if(k.equals("runs")){runs=v;}
			else if(k.equals("models")){models=v;}else if(k.equals("verified")){verified=v;}
			else if(k.equals("outprefix")){out=v;}else if(k.equals("expected")){expected=Integer.parseInt(v);}
			else if(k.equals("expectedmodels")){expectedModels=Integer.parseInt(v);}
			else if(k.equals("prokinputs")){prokInputs=Parse.parseBoolean(v);}else{throw new IllegalArgumentException(k);}
		}
		require(inputs!=null && sources!=null && contigs!=null && contigSources!=null && features!=null
			&& runs!=null && models!=null && verified!=null && out!=null && expected>0 && expectedModels>0,
			"Bind inputs, sources, contigs, contigsources, features, runs, models, verified and fresh outprefix");
		Shared.setThreads(1);ReadWrite.setZipThreads(1);ReadWrite.ZIPLEVEL=1;
		final EukPanelCmEvidence tool=new EukPanelCmEvidence();
		tool.prokInputs=prokInputs;
		tool.verified=EukPanelGffValidation.hashes(verified);
		for(String path:new String[]{inputs,sources,contigs,contigSources,features,runs,models}){tool.bound(path);}
		tool.panel.hashes=EukPanelGffValidation.hashes(sources);
		tool.panel.contigHashes=EukPanelGffValidation.hashes(contigSources);
		if(prokInputs){tool.loadProkInputs(inputs,expected);tool.loadProkContigs(contigs);}
		else{
			tool.panel.loadInputs(inputs,expected);tool.panel.loadContigs(contigs);
			read(inputs,"compartment\tassembly\tbin\ttaxid\tbases\tsequence_records\tfasta\tselected_fasta_sha80\tgff\tannotation_status\tgbff\tgbff_record_accession\tmetadata_status\tannotation_source_fasta\tpairing_provenance",p->{
				require(p.terms()==15,"Bound panel input schema changed");tool.sourceFastas.put(p.parseString(1),p.parseString(6));
			});
		}
		tool.loadModels(models);tool.loadRuns(runs,expectedModels);tool.write(features,out);
		Shared.closeStream(pp.outstream);
	}
	String bound(String path){
		final String hash=verified.get(path);require(hash!=null,"Input absent from producer-verified manifest: "+path);return hash;
	}
	void bound(String path,String hash){require(bound(path).equals(hash),"Declared and producer-verified bytes differ: "+path);}
	/** Loads the validator's native prok schema; domain is not an organelle compartment. */
	void loadProkInputs(String path,int expected){
		assert(prokInputs && panel.assemblies.isEmpty()):"Prok metadata must not be mixed with a preloaded euk panel";
		final ObjectSet<String> paths=new ObjectSet<String>(String.class);
		read(path,PROK_INPUT_HEADER,p->{
			require(p.terms()==12 && p.parseInt(0)==panel.assemblies.size()+1,"Prok input rows must retain the complete selected index order");
			final String id=p.parseString(1), domain=p.parseString(2), fasta=p.parseString(7), fastaHash=p.parseString(8), gff=p.parseString(9), gffHash=p.parseString(10);
			final int taxid=p.parseInt(3);final long selected=p.parseLong(4), observed=p.parseLong(5), records=p.parseLong(6);
			require(id.matches("GCF_[0-9]+\\.[0-9]+") && (domain.equals("bacteria") || domain.equals("archaea"))
				&& taxid>0 && selected>0 && selected==observed && records>0 && SeedOffsetTable.validSha80(p.parseString(11)),
				"Prok assembly requires explicit domain, positive observed dimensions and original download provenance");
			require(paths.add(fasta) && paths.add(gff),"Each assembly needs distinct FASTA/GFF paths");
			bound(fasta,fastaHash);bound(gff,gffHash);
			require(fastaHash.equals(panel.hashes.get(fasta)) && fastaHash.equals(panel.contigHashes.get(fasta))
				&& gffHash.equals(panel.hashes.get(gff)),"Prok input and contig observation must bind identical verified source bytes");
			final EukPanelGffValidation.Assembly a=new EukPanelGffValidation.Assembly(domain,id,gff,observed,records);a.gffHash=gffHash;
			require(panel.byAssembly.put(id,a)==null,"Duplicate prok assembly: "+id);
			panel.assemblies.add(a);sourceFastas.put(id,fasta);prokBindings.put(id,new ProkBinding(taxid,fasta,fastaHash));
		});
		require(panel.assemblies.size()==expected,"Prok input assembly count differs from declared panel size");
	}
	void loadProkContigs(String path){
		assert(prokInputs && !panel.assemblies.isEmpty()):"Observed contigs must belong to the selected prok assemblies";
		final int[] nativeTaxid=new int[1];
		read(path,PROK_CONTIG_HEADER,p->{
			require(p.terms()==10,"Prok contig observation schema changed");
			final EukPanelGffValidation.Assembly a=panel.byAssembly.get(p.parseString(1));
			require(a!=null && p.termEquals(a.group,0),"Observed prok contig domain/assembly differs from selected inputs");
			final ProkBinding binding=prokBindings.get(a.id);
			require(p.termEquals(binding.fasta,8) && p.termEquals(binding.hash,9) && p.parseInt(6)==binding.taxid,
				"Observed prok contig must bind its own assembly FASTA/hash and selected TaxID");
			final String id=p.parseString(2), bare=EukPanelCompartments.accession(id,nativeTaxid);final long length=p.parseLong(3);
			final String taxonomy=nativeTaxid[0]==0?"NOT_TAGGED":nativeTaxid[0]==binding.taxid?"MATCH":"MISMATCH";
			require(length>0 && p.termEquals("PROKARYOTIC",4) && p.parseInt(5)==nativeTaxid[0] && p.termEquals(taxonomy,7),
				"Prok contig length, native tag and taxonomy labels must remain observed facts, not nuclear labels");
			final EukPanelGffValidation.Contig c=new EukPanelGffValidation.Contig(id,length,"PROKARYOTIC",nativeTaxid[0]);
			require(a.literal.put(id,c)==null && !a.aliases.contains(bare),"Duplicate or ambiguous prok contig identity within "+a.id);
			EukPanelGffValidation.alias(a,bare,c,"VALIDATED_NATIVE_TAG_ACCESSION");a.observedRecords++;a.observedBases+=length;
		});
		for(EukPanelGffValidation.Assembly a:panel.assemblies){require(a.records==a.observedRecords && a.bases==a.observedBases,
			"Prok contig observations must conserve every selected assembly: "+a.id);}
	}
	void loadModels(String path){
		assert(models.size()==0):"One declared model roster binds all scan recipes";
		read(path,MODEL_HEADER,p->{
			require(p.terms()==8,"Model metadata requires eight fields");
			final Model m=new Model(p.parseString(0),p.parseString(1),p.parseString(2),p.parseString(3),
				p.parseInt(4),p.parseDouble(5),p.parseString(6),p.parseString(7));
			require(SeedOffsetTable.validSha80(m.bundle) && SeedOffsetTable.validSha80(m.sourceHash)
				&& m.accession.matches("RF[0-9]{5}") && m.clen>0 && Double.isFinite(m.ga),"Invalid pinned model metadata");
			bound(m.source,m.sourceHash);
			require(models.put(m.bundle+":"+m.accession,m)==null,"Duplicate CM-bytes/accession model key");
			modelOrder.add(m);
		});
		require(!modelOrder.isEmpty(),"An empty model roster cannot establish scan coverage");
	}
	void loadRuns(String path,int expectedModels){
		assert(runs.isEmpty()):"A scan is consumed once; controls/retries must be selected explicitly";
		final ObjectSet<String> rosterSet=new ObjectSet<String>(String.class);final ArrayList<String> roster=new ArrayList<String>();
		for(Model m:modelOrder){if(rosterSet.add(m.accession)){roster.add(m.accession);}}
		require(roster.size()==expectedModels,"Exact global accession roster differs from expected model count");
		final ObjectSet<String> ids=new ObjectSet<String>(String.class), tables=new ObjectSet<String>(String.class);
		read(path,RUN_HEADER,p->{
			require(p.terms()==16,"Run metadata requires sixteen fields");
			final Run r=new Run(p.parseString(0),p.parseString(1),p.parseString(2),p.parseString(3),
				p.parseString(4),p.parseString(5),p.parseString(6),p.parseString(7),p.parseString(8),
				p.parseString(9),p.parseString(10),p.parseString(11),p.parseString(12),p.parseString(13),p.parseString(14),p.parseString(15));
			require(ids.add(r.id) && tables.add(r.table),"Duplicate run ID or raw table would inflate evidence");
			final EukPanelGffValidation.Assembly a=panel.byAssembly.get(r.assembly);
			require(a!=null,"Scan assembly absent from selected panel: "+r.assembly);
			require(r.source.equals(sourceFastas.get(r.assembly)) && r.sourceHash.equals(panel.hashes.get(r.source)),
				"Run must bind this assembly's original FASTA, not just a target with matching taxid");
			require(r.format.equals("INFERNAL_1_1_5_TBLOUT17"),"Declare the actual raw format; derived22 is not tblout");
			bound(r.table,r.tableHash);bound(r.cm,r.cmHash);bound(r.scan,r.scanHash);
			bound(r.executable,r.executableHash);bound(r.command,r.commandHash);bound(r.seal);
			for(Model m:modelOrder){if(m.bundle.equals(r.cmHash)){
				final String key=r.assembly+":"+m.accession;
				require(scanned.put(key,m)==null,"More than one selected run covers assembly/accession: "+key);
				r.models.add(m);
			}}
			require(!r.models.isEmpty(),"Scan CM bundle has no declared metadata: "+r.id);
			runs.add(r);loadHits(r,a);
		});
		for(EukPanelGffValidation.Assembly a:panel.assemblies){
			for(String accession:roster){require(scanned.contains(a.id+":"+accession),"Assembly missing declared accession: "+a.id+" "+accession);}
			int n=0;for(Run r:runs){if(r.assembly.equals(a.id)){n+=r.models.size();}}
			require(n==expectedModels,"Every assembly needs the declared model coverage, including zero-hit scans: "+a.id+" "+n+"/"+expectedModels);
		}
		for(Bucket b:buckets){b.finish();}
	}
	void loadHits(Run run,EukPanelGffValidation.Assembly assembly){
		final ByteFile in=RrnaResourceIO.open(run.table);final RrnaCmTruthBuilder.Fields f=new RrnaCmTruthBuilder.Fields();
		boolean program=false,version=false,done=false;long lineNo=0;
		try{for(byte[] line;(line=in.nextLine())!=null;){
			lineNo++;if(line.length==0){continue;}require(!done,"Content after Infernal completion footer: "+run.id);
			if(line[0]=='#'){
				final String s=Euk18sShredCmSummary.text(line);
				if(s.startsWith("# Program:")){require(s.substring(10).trim().equals("cmsearch"),"Expected cmsearch");program=true;}
				if(s.startsWith("# Version:")){require(s.substring(10).trim().startsWith("1.1.5 "),"Expected pinned Infernal1.1.5");version=true;}
				if(s.equals("# [ok]")){done=true;}continue;
			}
			f.set(line);require(f.eq(4,"cm"),"An HMM-only row is not a CM hit");
			final Model m=models.get(run.cmHash+":"+f.text(3));require(m!=null,"Undeclared model in raw scan: "+run.id);
			final Hit h=new Hit(run,m,lineNo,f);
			require(h.mf>=1 && h.mt>=h.mf && h.mt<=m.clen && h.from>0 && h.to>0,"Raw model/sequence coordinate range invalid");
			require(h.strand.equals("+") && h.from<=h.to || h.strand.equals("-") && h.from>=h.to,"CM coordinate order must match strand");
			require(h.trunc.equals("no") || h.trunc.equals("5'") || h.trunc.equals("3'") || h.trunc.equals("5'&3'"),"Unknown truncation flag");
			require(h.inclusion.equals("!") || h.inclusion.equals("?"),"Unknown inclusion flag");
			require(Double.isFinite(h.score) && Double.isFinite(h.evalue) && h.evalue>=0,"Nonfinite CM score or invalid E-value");
			final EukPanelGffValidation.Resolution r=EukPanelGffValidation.resolve(assembly,h.target);
			h.mapping=r.kind;h.contig=r.contig;
			h.range=r.contig==null?"UNRESOLVED_SEQID":h.stop>r.contig.length?"OUT_OF_RANGE":"VALID_INTERVAL";
			run.hits.add(h);hits.add(h);
			if(h.range.equals("VALID_INTERVAL")){
				final String key=run.assembly+"\t"+h.contig.id+"\t"+h.strand;
				Bucket b=byContig.get(key);if(b==null){b=new Bucket();byContig.put(key,b);buckets.add(b);}b.hits.add(h);
			}
		}}finally{require(!in.close(),"Raw table read failed: "+run.table);}
		require(program && version && done,"Incomplete or unrecognized raw cmsearch table: "+run.table);
	}
	void write(String features,String prefix){
		assert(!runs.isEmpty()):"A completed scan with zero hits is still an explicit run";
		RrnaResourceIO.fresh(prefix,".hits.tsv.gz",".associations.tsv.gz",".sources.tsv.gz",".hit_overlaps.tsv.gz",".summary.tsv");
		final ArrayList<ByteStreamWriter> opened=new ArrayList<ByteStreamWriter>(5);
		final ByteBuilder b=new ByteBuilder();
		try{
			final ByteStreamWriter raw=openWriter(prefix+".hits.tsv.gz",opened), edges=openWriter(prefix+".associations.tsv.gz",opened);
			final ByteStreamWriter src=openWriter(prefix+".sources.tsv.gz",opened), cross=openWriter(prefix+".hit_overlaps.tsv.gz",opened), summary=openWriter(prefix+".summary.tsv",opened);
			raw.println("assembly\trun_id\ttable\ttable_sha80\ttable_line\tcm_sha80\tmodel\tfamily\tsubtype\ttarget\tnative_seqid\tcompartment\tmapping_status\trange_status\tmodel_from\tmodel_to\tclen\tmodel_span\tmodel_coverage\tseq_from\tseq_to\tstrand\ttrunc\tinclusion\tbits\tevalue\tga\tga_margin\thistorical_overlaps");
			edges.println("assembly\tsource_gff\tsource_line\tstatus\trun_id\ttable_line\tmodel\toverlap_bp\tsource_fraction\thit_fraction");
			src.println(FEATURE_HEADER+"\tassociation_status\toverlapping_hits\tdistinct_models");
			cross.println("assembly\tnative_seqid\tstrand\trun_a\tline_a\tmodel_a\trun_b\tline_b\tmodel_b\toverlap_bp\trelation");
			summary.println("assembly\trun_id\tmodel\tcm_sha80\traw_hits\tvalid_interval_hits\tunresolved_hits\tout_of_range_hits\tincluded_hits\tfull_model_span_hits\tuntruncated_hits\tga_or_higher_hits\thits_with_historical_overlap");
			associate(features,src,edges);
			for(Hit h:hits){
				b.clear().append(h.run.assembly).tab().append(h.run.id).tab().append(h.run.table).tab().append(h.run.tableHash).tab().append(h.line).tab().append(h.model.bundle).tab().append(h.model.accession).tab().append(h.model.family).tab().append(h.model.subtype).tab().append(h.target).tab().append(h.contig==null?"NA":h.contig.id).tab().append(h.contig==null?"UNKNOWN":h.contig.compartment).tab().append(h.mapping).tab().append(h.range);
				b.tab().append(h.mf).tab().append(h.mt).tab().append(h.model.clen).tab().append(h.mt-h.mf+1).tab().append((h.mt-h.mf+1)/(double)h.model.clen,9).tab().append(h.from).tab().append(h.to).tab().append(h.strand).tab().append(h.trunc).tab().append(h.inclusion).tab().append(h.score,9).tab().append(h.evalueText).tab().append(h.model.ga,9).tab().append(h.score-h.model.ga,9).tab().append(h.associations);raw.print(b.nl());
			}
			for(Bucket bucket:buckets){for(int i=0;i<bucket.hits.size();i++){
				final Hit a=bucket.hits.get(i);
				for(int j=i+1;j<bucket.hits.size();j++){
					final Hit z=bucket.hits.get(j);if(z.start>a.stop){break;}
					final long overlap=overlap(a.start,a.stop,z.start,z.stop);
					b.clear().append(a.run.assembly).tab().append(a.contig.id).tab().append(a.strand).tab().append(a.run.id).tab().append(a.line).tab().append(a.model.accession).tab().append(z.run.id).tab().append(z.line).tab().append(z.model.accession).tab().append(overlap).tab().append(a.model.accession.equals(z.model.accession)?"SAME_MODEL":"CROSS_MODEL");cross.print(b.nl());
				}
			}}
			for(Run run:runs){for(Model m:run.models){
				final long[] c=new long[9];for(Hit h:run.hits){if(h.model!=m){continue;}
					c[0]++;c[h.range.equals("VALID_INTERVAL")?1:h.range.equals("UNRESOLVED_SEQID")?2:3]++;
					if(h.inclusion.equals("!")){c[4]++;}if(h.mf==1 && h.mt==m.clen){c[5]++;}
					if(h.trunc.equals("no")){c[6]++;}if(h.score>=m.ga){c[7]++;}if(h.associations>0){c[8]++;}
				}
				require(c[0]==c[1]+c[2]+c[3],"Every raw hit has exactly one interval status");
				b.clear().append(run.assembly).tab().append(run.id).tab().append(m.accession).tab().append(m.bundle);for(long n:c){b.tab().append(n);}summary.print(b.nl());
			}}
		}finally{RrnaResourceIO.finish(opened.toArray(new ByteStreamWriter[opened.size()]));}
		System.err.println("CM_EVIDENCE_COMPLETE assemblies="+panel.assemblies.size()+" runs="+runs.size()+" raw_hits="+hits.size()+" historical_rows="+sourceRows+" associations="+associationRows+" no_truth_admission=true");
	}
	void associate(String path,ByteStreamWriter sources,ByteStreamWriter edges){
		final ByteBuilder b=new ByteBuilder();final ObjectMap<String,SourceState> states=new ObjectMap<String,SourceState>(String.class,SourceState.class);
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		SourceState current=null;String lastContig=null;Bucket plus=null,minus=null;
		try{RrnaResourceIO.header(in,FEATURE_HEADER,path);for(byte[] line;(line=in.nextLine())!=null;){
			p.set(line);require(p.terms()==23,"Historical sidecar must preserve all nine source fields");
			if(current==null || !p.termEquals(current.a.id,1)){
				final String key=p.parseString(1);current=states.get(key);
				if(current==null){final EukPanelGffValidation.Assembly a=panel.byAssembly.get(key);require(a!=null,"Historical assembly absent from panel");current=new SourceState(a);states.put(key,current);}
				lastContig=null;plus=null;minus=null;
			}
			final String assembly=current.a.id,gff=current.a.gff;final long sourceLine=p.parseLong(4);
			require(p.termEquals(gff,2) && p.termEquals(current.a.gffHash,3),"Historical row differs from bound source identity");
			if(prokInputs){checkProkFeature(current.a,p);}
			require(sourceLine>current.lastLine,"Historical source lines must be unique and ordered within assembly");current.lastLine=sourceLine;
			final boolean valid=p.termEquals("VALID_INTERVAL",5), forward=p.termEquals("+",20), reverse=p.termEquals("-",20);
			long count=0;int distinct=0;final long epoch=sourceRows+1;
			String outcome;
			if(valid){outcome=forward || reverse?"NO_SAME_STRAND_OVERLAP":"STRAND_UNSPECIFIED";}
			else if(p.termEquals("OUT_OF_RANGE",5)){outcome="NOT_ASSOCIATED_OUT_OF_RANGE";}
			else if(p.termEquals("UNRESOLVED_SEQID",5)){outcome="NOT_ASSOCIATED_UNRESOLVED_SEQID";}
			else{require(p.termEquals("INVALID_COORDINATES",5),"Unknown source-validation status");outcome="NOT_ASSOCIATED_INVALID_COORDINATES";}
			if(valid && (forward || reverse)){
				final long start=p.parseLong(17), stop=p.parseLong(18);require(start>0 && stop>=start && stop<=p.parseLong(9),"Valid sidecar interval contradicts coordinates");
				if(lastContig==null || !p.termEquals(lastContig,8)){
					lastContig=p.parseString(8);plus=byContig.get(assembly+"\t"+lastContig+"\t+");minus=byContig.get(assembly+"\t"+lastContig+"\t-");
				}
				final Bucket bucket=forward?plus:minus;
				if(bucket!=null){for(int i=bucket.firstEndingAtLeast(start);i<bucket.hits.size();i++){
					final Hit h=bucket.hits.get(i);if(h.start>stop){break;}final long n=overlap(start,stop,h.start,h.stop);if(n==0){continue;}
					count++;h.associations++;associationRows++;if(h.model.lastSource!=epoch){h.model.lastSource=epoch;distinct++;}
					b.clear().append(assembly).tab().append(gff).tab().append(sourceLine).append("\tOVERLAP\t").append(h.run.id).tab().append(h.line).tab().append(h.model.accession).tab().append(n).tab().append(n/(double)(stop-start+1),9).tab().append(n/(double)(h.stop-h.start+1),9);edges.print(b.nl());
				}}
			}
			if(count>0){outcome=distinct>1?"MULTIPLE_MODELS":count>1?"MULTIPLE_HITS":"ONE_HIT";}
			else{b.clear().append(assembly).tab().append(gff).tab().append(sourceLine).tab().append(outcome).append("\tNA\tNA\tNA\tNA\tNA\tNA");edges.print(b.nl());}
			sources.print(b.clear().append(line).tab().append(outcome).tab().append(count).tab().append(distinct).nl());sourceRows++;
		}}finally{require(!in.close(),"Historical sidecar read failed");}
	}
	void checkProkFeature(EukPanelGffValidation.Assembly a,LineParser1 p){
		assert(prokInputs):"Only the explicit prok schema uses this additional observation binding";
		require(p.termEquals(a.group,0),"Historical prok feature domain differs from bound assembly");
		final EukPanelGffValidation.Resolution r=EukPanelGffValidation.resolve(a,p.parseString(14));
		require(p.termEquals(r.kind,6) && p.parseInt(11)==r.taxid,"Historical prok feature mapping/native tag differs from its original seqid");
		if(r.contig==null){require(p.termEquals("NA",8) && p.termEquals("NA",9) && p.termEquals("UNKNOWN",10),"Unresolved prok feature must not invent a contig");}
		else{require(p.termEquals(r.contig.id,8) && p.parseLong(9)==r.contig.length && p.termEquals("PROKARYOTIC",10)
			&& p.parseInt(12)==r.contig.taxid,"Historical prok feature contradicts observed contig identity/length/compartment");}
	}
	static ByteStreamWriter openWriter(String path,ArrayList<ByteStreamWriter> opened){
		assert(path!=null && opened!=null):"Every started writer must be closed even if a later output fails to open";
		final ByteStreamWriter w=RrnaResourceIO.writer(path);opened.add(w);return w;
	}
	static long overlap(long a,long z,long b,long y){return Math.max(0,Math.min(z,y)-Math.max(a,b)+1);}
	static void read(String path,String header,Rows rows){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{RrnaResourceIO.header(in,header,path);for(byte[] line;(line=in.nextLine())!=null;){p.set(line);rows.row(p);}}
		finally{require(!in.close(),"Metadata read failed: "+path);}
	}
	interface Rows{void row(LineParser1 p);}
	public static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	public static final String MODEL_HEADER="bundle_sha80\taccession\tfamily\tsubtype\tclen\tga\tsource_cm\tsource_cm_sha80";
	public static final String RUN_HEADER="run_id\tassembly\traw_table\traw_sha80\tcm_file\tcm_sha80\tscan_fasta\tscan_sha80\texecutable\texecutable_sha80\tcommand\tcommand_sha80\tresult_manifest\tsource_fasta\tsource_sha80\tformat";
	public static final String FEATURE_HEADER="selection_group\tassembly\tsource_gff\tsource_sha80\tsource_line\tvalidation_status\tmapping_status\trange_status\tnative_seqid\tcontig_length\tcompartment\tgff_taxid\tnative_taxid\ttaxonomy_status\toriginal_seqid\toriginal_source\toriginal_type\toriginal_start\toriginal_end\toriginal_score\toriginal_strand\toriginal_phase\toriginal_attributes";
	public static final String PROK_INPUT_HEADER="index\tassembly\tdomain\tselected_taxid\tselected_bases\tobserved_bases\tcontigs\tfasta\tfasta_sha80\tgff\tgff_sha80\tdownload_seal";
	public static final String PROK_CONTIG_HEADER="domain\tassembly\tseqid\tlength\tcompartment\tnative_taxid\tselected_taxid\ttaxonomy_status\tfasta\tfasta_sha80";
	static final class ProkBinding{ProkBinding(int t,String f,String h){taxid=t;fasta=f;hash=h;}final int taxid;final String fasta,hash;}
	static final class Model{
		Model(String b,String a,String f,String s,int n,double g,String p,String h){bundle=b;accession=a;family=f;subtype=s;clen=n;ga=g;source=p;sourceHash=h;}
		final String bundle,accession,family,subtype,source,sourceHash;final int clen;final double ga;long lastSource;
	}
	static final class Run{
		Run(String i,String a,String t,String th,String c,String ch,String s,String sh,String e,String eh,String p,String ph,String seal_,String source_,String sourceHash_,String format_){id=i;assembly=a;table=t;tableHash=th;cm=c;cmHash=ch;scan=s;scanHash=sh;executable=e;executableHash=eh;command=p;commandHash=ph;seal=seal_;source=source_;sourceHash=sourceHash_;format=format_;}
		final String id,assembly,table,tableHash,cm,cmHash,scan,scanHash,executable,executableHash,command,commandHash,seal,source,sourceHash,format;
		final ArrayList<Model> models=new ArrayList<Model>();final ArrayList<Hit> hits=new ArrayList<Hit>();
	}
	static final class Hit implements Comparable<Hit>{
		Hit(Run r,Model m,long l,RrnaCmTruthBuilder.Fields f){
			run=r;model=m;line=l;target=f.text(0);mf=f.number(5);mt=f.number(6);from=f.number(7);to=f.number(8);
			start=Math.min(from,to);stop=Math.max(from,to);strand=f.text(9);trunc=f.text(10);inclusion=f.text(16);
			score=Parse.parseDouble(f.bytes,f.start[14],f.end[14]);evalueText=f.text(15);evalue=Double.parseDouble(evalueText);
		}
		@Override public int compareTo(Hit other){int c=Integer.compare(start,other.start);if(c==0){c=Integer.compare(stop,other.stop);}if(c==0){c=run.id.compareTo(other.run.id);}return c==0?Long.compare(line,other.line):c;}
		final Run run;final Model model;final long line;final String target,strand,trunc,inclusion,evalueText;final int mf,mt,from,to,start,stop;final double score,evalue;
		String mapping,range;EukPanelGffValidation.Contig contig;long associations;
	}
	static final class Bucket{
		void finish(){Collections.sort(hits);maxEnd=new int[hits.size()];int max=0;for(int i=0;i<hits.size();i++){max=Math.max(max,hits.get(i).stop);maxEnd[i]=max;}}
		int firstEndingAtLeast(long start){assert(maxEnd!=null):"Prefix endpoint index is built after sorting retained hits";int lo=0,hi=maxEnd.length;while(lo<hi){final int mid=(lo+hi)>>>1;if(maxEnd[mid]<start){lo=mid+1;}else{hi=mid;}}return lo;}
		final ArrayList<Hit> hits=new ArrayList<Hit>();int[] maxEnd;
	}
	static final class SourceState{SourceState(EukPanelGffValidation.Assembly a_){a=a_;}final EukPanelGffValidation.Assembly a;long lastLine;}
	final EukPanelGffValidation panel=new EukPanelGffValidation();
	ObjectMap<String,String> verified;
	boolean prokInputs;
	final ObjectMap<String,ProkBinding> prokBindings=new ObjectMap<String,ProkBinding>(String.class,ProkBinding.class);
	final ObjectMap<String,String> sourceFastas=new ObjectMap<String,String>(String.class,String.class);
	final ObjectMap<String,Model> models=new ObjectMap<String,Model>(String.class,Model.class), scanned=new ObjectMap<String,Model>(String.class,Model.class);
	final ObjectMap<String,Bucket> byContig=new ObjectMap<String,Bucket>(String.class,Bucket.class);
	final ArrayList<Model> modelOrder=new ArrayList<Model>();final ArrayList<Run> runs=new ArrayList<Run>();
	final ArrayList<Hit> hits=new ArrayList<Hit>();final ArrayList<Bucket> buckets=new ArrayList<Bucket>();long sourceRows,associationRows;
}
