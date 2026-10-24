package prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.locks.ReadWriteLock;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import idaligner.AlignmentStats;
import idaligner.QuantumAligner;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import structures.ByteBuilder;
import structures.IntList;
import structures.IntHashMap;
import tax.TaxTree;
import template.Accumulator;
import template.ThreadWaiter;

/** Coverage-selected mitochondrial candidate models, never production promotion.
 * Every original copy survives exact-sequence caching. Alignment uses the shorter
 * sequence as complete query; this differs from complete-model caller verification.
 * @author Brian, Ganyu */
public final class NcrnaMitoCoverageSelect implements Accumulator<NcrnaMitoCoverageSelect.Worker> {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args,NcrnaMitoCoverageSelect.class,false);
		try{final NcrnaMitoCoverageSelect x=new NcrnaMitoCoverageSelect();x.parse(pp.args);x.run();}
		finally{Shared.closeStream(pp.outstream);}
	}
	private void parse(String[] args){
		for(String arg:args){final String[] kv=Parse.splitOnFirst(arg,'=');require(kv.length==2,"Use flag=value");
			final String key=kv[0].toLowerCase(Locale.ROOT);require(KNOWN.contains("|"+key+"|") && flags.put(key,kv[1])==null,"Unknown/repeated flag: "+key);}
		for(String key:INPUTS){require(flags.get(key)!=null,"Missing bound input: "+key);}
		prefix=flags.get("outprefix");mode=value("mode","pilot");require(prefix!=null && (mode.equals("pilot") || mode.equals("full") || mode.equals("select")),"Use explicit mode=pilot|full|select and a fresh outprefix");
		corpusExpected=integer("expected",16110);panelExpected=integer("expectedpanel",22);ssuExpected=integer("expectedssu",155);nohitExpected=integer("expectednohit",142);
		pilotSize=integer("pilotsize",64);threads=integer("t",1);cut=Float.parseFloat(value("cut","0.70"));
		require(corpusExpected>0 && panelExpected>0 && ssuExpected>0 && nohitExpected>0 && pilotSize>0 && threads>0 && Float.isFinite(cut) && cut>0 && cut<=1,"Positive denominators/threads and an explicit valid coverage threshold required");
		pins=EukPanelGffValidation.hashes(required("verified"));for(String key:INPUTS){require(pins.get(flags.get(key))!=null,"Producer must verify input "+key+": "+flags.get(key));}
		if(mode.equals("select")){require(flags.get("matrix")!=null && pins.get(flags.get("matrix"))!=null,"Selection requires a producer-verified complete matrix");}
		Shared.setThreads(threads);ReadWrite.setZipThreads(Math.min(4,threads));Shared.TRIM_READ_DESCRIPTION=false;Read.TO_UPPER_CASE=true;
		RrnaResourceIO.fresh(prefix,".matrix.tsv.gz",".unique.fa.gz",".members.tsv.gz",".models.tsv",".selection.tsv",".selected.fa",".composition.tsv",".summary.tsv",".strata.tsv",".pilot.tsv");
	}
	private void run(){
		final boolean old=TaxTree.LOAD_FAKE_TREE;try{TaxTree.LOAD_FAKE_TREE=false;tree=TaxTree.loadTaxTreeFile(required("tree"));}finally{TaxTree.LOAD_FAKE_TREE=old;}
		loadModels("ssu",required("ssumodels"),ssuExpected);loadModels("nohit",required("nohitmodels"),nohitExpected);
		loadSourceMembers("ssu",required("ssumembers"));loadSourceMembers("nohit",required("nohitmembers"));
		loadCorpus();loadPanel();require(members.size()==corpusExpected+panelExpected,"Every construction copy and panel locus must survive loading");
		chooseMeasured();
		for(Sequence s:unique){if(s.measured){measuredUnique++;}}
		// Freeze the workload and its sequence/model bindings before any worker starts.
		writeBindings();if(mode.equals("pilot")){writePilot();}
		if(mode.equals("select")){loadMatrix(required("matrix"));}
		else{
			final ArrayList<Worker> workers=new ArrayList<Worker>();for(int i=0;i<Math.min(threads,unique.size());i++){workers.add(new Worker(this));}
			require(ThreadWaiter.startAndWait(workers,this),"An alignment worker failed; incomplete coverage cannot select a library");
			require(alignmentCalls==(long)measuredUnique*canonicalModels,"Completed workers must align each measured unique sequence against each canonical model exactly once");writeMatrix();
		}
		if(mode.equals("pilot")){
			System.err.println("MITO_COVERAGE_PILOT_PASS copies="+members.size()+" all_unique="+unique.size()+" measured_unique="+measuredUnique+" models="+models.size()+" canonical_models="+canonicalModels+" actual_alignment_calls="+alignmentCalls+" matrix_rows="+((long)measuredUnique*models.size())+" no_resource_selection=true");return;
		}
		select();
	}
	private void loadModels(String role,String path,int expected){
		final int before=models.size();
		RrnaResourceIO.stream(path,r->{final String token=RrnaResourceIO.first(r.id),key=role+":"+token;
			require(r.bases!=null && r.length()>0 && modelByKey.get(key)==null && token.matches("[A-Za-z0-9_]+"),"Unique safe source model ID and sequence required: "+key);
			final Model m=new Model(key,role,token,r.id,r.bases.clone(),models.size());
			m.canonical=m.index;for(Model earlier:models){if(Arrays.equals(earlier.bases,m.bases)){m.canonical=earlier.canonical;break;}}
			if(m.canonical==m.index){canonicalModels++;}
			models.add(m);modelByKey.put(key,m);
		});require(models.size()-before==expected,"Complete frozen model count changed for "+role);
	}
	private void loadSourceMembers(String role,String path){
		final ByteFile input=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');int count=0;
		try{p.set(input.nextLine());require(p.terms()==22 && p.termEquals("locus_id",1) && p.termEquals("assigned_model",16),"Expected native cluster-evaluation membership schema");
			for(byte[] line;(line=input.nextLine())!=null;){count++;p.set(line);require(p.terms()==22 && p.termEquals(role,0),"Frozen cluster membership role/schema changed");
				final Source s=new Source(p.parseString(1),p.parseString(2),p.parseInt(3),p.parseInt(4),bool(p.parseString(5)),bool(p.parseString(6)),role,path);
				require(sourceMembers.put(s.id,s)==null,"Construction source locus must be globally unique across source cohorts: "+s.id);
				if(!p.termEquals("NA",16)){
					final Model m=modelByKey.get(role+":"+RrnaResourceIO.first(p.parseString(16)));
					require(m!=null,"Assigned source cluster absent from frozen candidate pool: "+p.parseString(16));s.model=m.key;
					composition(m.sourceComposition,phylum(s.taxid),s.uid);
				}
			}
		}finally{require(!input.close(),"Source cluster membership read failed");}
		require(count==integer(role.equals("ssu")?"expectedssucopies":"expectednohitcopies",role.equals("ssu")?11420:7050),"Full source-cluster membership denominator changed: "+role);
	}
	private void loadCorpus(){
		final ObjectSet<String> seen=new ObjectSet<String>(String.class);
		RrnaResourceIO.stream(required("corpus"),r->{final String id=RrnaResourceIO.first(r.id);final Source source=sourceMembers.get(id);
			require(source!=null && seen.add(id) && r.length()==source.length && TaxTree.parseTaxID(r.id)==source.taxid
				&& RrnaResourceIO.field(r.id,"unique").equals(source.uid) && bool(RrnaResourceIO.field(r.id,"partial"))==source.partial
				&& bool(RrnaResourceIO.field(r.id,"pseudo"))==source.pseudo,"Corpus must exactly bind native source membership: "+id);
			final Sequence s=sequence(r.bases);s.construction=true;
			final Sequence prior=sourceSequences.get(source.uid);require(prior==null || prior==s,"One original unique ID cannot bind conflicting sequences: "+source.uid);sourceSequences.put(source.uid,s);
			members.add(new Member("CONSTRUCTION",id,s,source.taxid,source.partial,source.pseudo,source.model,source.role,source.path));
		});require(members.size()==corpusExpected,"All original construction copies are required");
	}
	private void loadPanel(){
		final ObjectMap<String,Panel> panel=new ObjectMap<String,Panel>(String.class,Panel.class),byTruth=new ObjectMap<String,Panel>(String.class,Panel.class);
		ByteFile input=RrnaResourceIO.open(required("panelcases"));final LineParser1 p=new LineParser1('\t');
		try{p.set(input.nextLine());require(p.terms()>=40 && p.termEquals("host_taxid",37) && p.termEquals("audit_id",29),"Expected source-bound case table with host taxonomy");
			for(byte[] line;(line=input.nextLine())!=null;){p.set(line);require(p.terms()>=40,"Truncated panel case");if(!p.termEquals("ORGANELLE_SSU_TRUTH",1)){continue;}
				final Panel a=new Panel(p);require(panel.put(a.id,a)==null && byTruth.put(a.truthId,a)==null,"Panel case and frozen truth IDs must be unique");}
		}finally{require(!input.close(),"Panel cases read failed");}
		require(panel.size()==panelExpected,"Complete panel locus population required");input=RrnaResourceIO.open(required("truth"));
		try{RrnaResourceIO.header(input,EukPanelTruthFreeze.HEADER,required("truth"));
			for(byte[] line;(line=input.nextLine())!=null;){p.set(line);require(p.terms()==35,"Frozen truth schema changed");final Panel a=byTruth.get(p.parseString(0));if(a==null){continue;}
				require(!a.bound && p.termEquals(a.assembly,1) && p.termEquals(a.contig,2) && p.termEquals(a.strand,3)
					&& p.parseInt(4)==a.start && p.parseInt(5)==a.stop && p.termEquals("SSU",6),"Panel span must bind its exact frozen SSU truth: "+a.id);
				a.partial=!bit(p,15) || bit(p,16);a.geometry=p.parseString(27);a.bound=true;}
		}finally{require(!input.close(),"Frozen panel truth read failed");}
		RrnaResourceIO.stream(required("panelgenes"),r->{final Panel a=panel.get(RrnaResourceIO.first(r.id));if(a==null){return;}
			require(a.bound && !a.seen && r.length()==a.stop-a.start+1 && RrnaResourceIO.field(r.id,"assembly").equals(a.assembly)
				&& RrnaResourceIO.field(r.id,"contig").equals(a.contig) && RrnaResourceIO.field(r.id,"start").equals(Integer.toString(a.start))
				&& RrnaResourceIO.field(r.id,"stop").equals(Integer.toString(a.stop)) && RrnaResourceIO.field(r.id,"strand").equals(a.strand),"Panel exact-gene sequence must match exported source geometry: "+a.id);
			a.seen=true;final Sequence s=sequence(r.bases);s.panel=true;members.add(new Member("PANEL_INCLUDED",a.id,s,a.taxid,a.partial,false,"NA",a.geometry,required("truth")));});
		for(Panel a:panel.values()){if(a!=null){require(a.seen,"Missing bound panel gene: "+a.id);}}
	}
	private Sequence sequence(byte[] bases){
		require(bases!=null && bases.length>0,"Exact caching cannot admit an empty sequence");final Key key=new Key(bases);Sequence s=cache.get(key);
		if(s==null){s=new Sequence(bases.clone(),unique.size(),models.size());cache.put(new Key(s.bases),s);unique.add(s);}return s;
	}
	private void chooseMeasured(){
		if(!mode.equals("pilot")){for(Sequence s:unique){s.measured=true;}return;}
		final ArrayList<Sequence> ordered=new ArrayList<Sequence>();for(Sequence s:unique){if(s.construction){ordered.add(s);}if(s.panel){s.measured=true;}}
		ordered.sort((a,b)->{int n=Integer.compare(a.bases.length,b.bases.length);return n==0?Integer.compare(a.index,b.index):n;});
		final int count=Math.min(pilotSize,ordered.size());for(int i=0;i<count;i++){ordered.get(count==1?0:(int)((long)i*(ordered.size()-1)/(count-1))).measured=true;}
	}
	private void writeMatrix(){
		final ByteStreamWriter out=RrnaResourceIO.writer(prefix+".matrix.tsv.gz");final ByteBuilder b=new ByteBuilder();
		try{out.println("#MITO_COVERAGE_MATRIX_V1\tSHORTER_QUERY_QUANTUM\t"+mode);
			for(String key:INPUTS){out.print(b.clear().append("#input\t").append(key).tab().append(pins.get(flags.get(key))).nl());}
			out.println("unique_id\tmodel_id\tmember_length\tmodel_length\tquery_side\tidentity\ttrace_identity\tpaired_columns");
			for(Sequence s:unique){if(s.measured){for(Model m:models){out.print(b.clear().append(s.id).tab().append(m.key).tab().append(s.bases.length).tab().append(m.bases.length)
				.tab().append(s.bases.length<=m.bases.length?"MEMBER":"MODEL").tab().appendSlow(s.identity[m.index]).tab().appendSlow(s.trace[m.index]).tab().append(s.paired[m.index]).nl());}}}
		}finally{RrnaResourceIO.finish(out);}
	}
	private void loadMatrix(String path){
		final ByteFile input=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');long rows=0;
		try{RrnaResourceIO.header(input,"#MITO_COVERAGE_MATRIX_V1\tSHORTER_QUERY_QUANTUM\tfull",path);
			for(String key:INPUTS){RrnaResourceIO.header(input,"#input\t"+key+'\t'+pins.get(flags.get(key)),path);}
			RrnaResourceIO.header(input,"unique_id\tmodel_id\tmember_length\tmodel_length\tquery_side\tidentity\ttrace_identity\tpaired_columns",path);
			for(Sequence s:unique){for(Model m:models){final byte[] line=input.nextLine();require(line!=null,"Incomplete coverage matrix");p.set(line);
				require(p.terms()==8 && p.termEquals(s.id,0) && p.termEquals(m.key,1) && p.parseInt(2)==s.bases.length && p.parseInt(3)==m.bases.length
					&& p.termEquals(s.bases.length<=m.bases.length?"MEMBER":"MODEL",4),"Matrix source order/identity/length mismatch");
				s.identity[m.index]=p.parseFloat(5);s.trace[m.index]=p.parseFloat(6);s.paired[m.index]=p.parseInt(7);validate(s,m);rows++;}}
			require(input.nextLine()==null && rows==(long)unique.size()*models.size(),"Matrix must contain exactly every expected pair");
		}finally{require(!input.close(),"Coverage matrix read failed");}
	}
	private void writeBindings(){
		final ByteStreamWriter seq=RrnaResourceIO.writer(prefix+".unique.fa.gz"),out=RrnaResourceIO.writer(prefix+".models.tsv");final ByteBuilder b=new ByteBuilder();
		try{for(Sequence s:unique){seq.print(b.clear().append('>').append(s.id).append(" construction=").append(s.construction).append(" panel=").append(s.panel).append(" measured=").append(s.measured).nl().append(s.bases).nl());}
			out.println("model_id\toutput_id\tlength\tcanonical_model\tsource_namespace\toriginal_header");
			for(Model m:models){out.print(b.clear().append(m.key).tab().append(m.output).tab().append(m.bases.length).tab().append(models.get(m.canonical).key).tab().append(m.role).tab().append(m.header).nl());}
		}finally{RrnaResourceIO.finish(seq,out);}
	}
	private void writePilot(){
		final ByteStreamWriter out=RrnaResourceIO.writer(prefix+".pilot.tsv");final ByteBuilder b=new ByteBuilder();
		try{out.println("unique_id\tlength\tconstruction_member\tpanel_member\tplanned_alignments");for(Sequence s:unique){if(s.measured){out.print(b.clear().append(s.id).tab().append(s.bases.length).tab().append(s.construction).tab().append(s.panel).tab().append(canonicalModels).nl());}}}
		finally{RrnaResourceIO.finish(out);}
	}
	private void select(){
		final int words=(unique.size()+63)/64;final long[][] cover=new long[models.size()][words];
		for(Model m:models){if(m.canonical==m.index){for(Sequence s:unique){if(s.identity[m.index]>=cut){cover[m.index][s.index>>>6]|=1L<<(s.index&63);}}}}
		final int[] selected=greedy(cover,unique.size());final int[] assigned=new int[unique.size()];Arrays.fill(assigned,-1);
		final ByteStreamWriter steps=RrnaResourceIO.writer(prefix+".selection.tsv"),fasta=RrnaResourceIO.writer(prefix+".selected.fa"),rows=RrnaResourceIO.writer(prefix+".members.tsv.gz"),composition=RrnaResourceIO.writer(prefix+".composition.tsv"),summary=RrnaResourceIO.writer(prefix+".summary.tsv"),strata=RrnaResourceIO.writer(prefix+".strata.tsv");
		final ObjectMap<String,Composition> totals=new ObjectMap<String,Composition>(String.class,Composition.class),coveredGroups=new ObjectMap<String,Composition>(String.class,Composition.class);
		final ByteBuilder b=new ByteBuilder();int covered=0;long coveredCopies=0;int partial=0,panel=0;
		try{steps.println("step\tmodel_id\tnew_unique\tnew_copies\tcovered_unique\tremaining_unique");composition.println("model_id\tscope\tphylum\tcopies\tunique_sequences");
			for(int step=0;step<selected.length;step++){final int index=selected[step];final Model m=models.get(index);int gain=0,copies=0;
				for(Sequence s:unique){if(assigned[s.index]<0 && s.identity[index]>=cut){assigned[s.index]=index;gain++;}}
				final ObjectMap<String,Composition> groups=new ObjectMap<String,Composition>(String.class,Composition.class);
				for(Member member:members){if(assigned[member.sequence.index]==index){copies++;composition(groups,phylum(member.taxid),member.sequence.id);}}
				covered+=gain;steps.print(b.clear().append(step+1).tab().append(m.key).tab().append(gain).tab().append(copies).tab().append(covered).tab().append(unique.size()-covered).nl());
				fasta.print(b.clear().append('>').append(m.output).append(" source=").append(m.key).nl().append(m.bases).nl());
				for(Model alias:models){if(alias.canonical==index){writeComposition(composition,b,alias,"SOURCE_CLUSTER_FOR_SELECTED_"+m.key,alias.sourceComposition);}}
				writeComposition(composition,b,m,"NEWLY_COVERED_TARGETS",groups);
			}
			rows.println("role\tlocus_id\tunique_id\ttaxid\tphylum\tpartial\tpseudo\tstratum\tsource_model\tassigned_cover_model\tbest_selected_model\tselected_identity\tselected_trace_identity\tmember_paired_fraction\tmodel_paired_fraction\tbest_pool_model\tbest_pool_identity\tstatus\tprovenance");
			for(Member member:members){final Sequence s=member.sequence;int best=-1,pool=0;for(Model m:models){if(s.identity[m.index]>s.identity[pool]){pool=m.index;}}for(int i:selected){if(best<0 || s.identity[i]>s.identity[best]){best=i;}}
				final boolean yes=assigned[s.index]>=0;if(yes){coveredCopies++;}partial+=member.partial?1:0;panel+=member.role.equals("PANEL_INCLUDED")?1:0;
				final String stratum=member.role+":"+(member.partial?"PARTIAL":"NOT_MARKED_PARTIAL");composition(totals,stratum,s.id);if(yes){composition(coveredGroups,stratum,s.id);}
				b.clear().append(member.role).tab().append(member.id).tab().append(s.id).tab().append(member.taxid).tab().append(phylum(member.taxid)).tab().append(member.partial).tab().append(member.role.equals("PANEL_INCLUDED")?"NA":Boolean.toString(member.pseudo))
					.tab().append(member.stratum).tab().append(member.sourceModel).tab().append(yes?models.get(assigned[s.index]).key:"NA").tab().append(best<0?"NA":models.get(best).key).tab();
				if(best<0){b.append("NA\tNA\tNA\tNA");}else{b.appendSlow(s.identity[best]).tab().appendSlow(s.trace[best]).tab().append(s.paired[best]/(double)s.bases.length,9).tab().append(s.paired[best]/(double)models.get(best).bases.length,9);}
				rows.print(b.tab().append(models.get(pool).key).tab().appendSlow(s.identity[pool]).tab().append(yes?"COVERED_SHORTER_QUERY":"UNCOVERED_AT_CUT").tab().append(member.source).nl());
			}
			summary.println("copies\tunique_sequences\tmodels\tselected_models\tcut\tcovered_unique\tuncovered_unique\tcovered_copies\tuncovered_copies\tpartial_copies\tpanel_included_copies\tinterpretation\tmode\tcanonical_models\tactual_alignment_calls\tmatrix_rows");
			summary.print(b.clear().append(members.size()).tab().append(unique.size()).tab().append(models.size()).tab().append(selected.length).tab().appendSlow(cut).tab().append(covered).tab().append(unique.size()-covered)
				.tab().append(coveredCopies).tab().append(members.size()-coveredCopies).tab().append(partial).tab().append(panel).append("\tSHORTER_QUERY_COVERAGE_NOT_FULL_MODEL_CALLER_VERIFICATION")
				.tab().append(mode).tab().append(canonicalModels).tab().append(alignmentCalls).tab().append((long)measuredUnique*models.size()).nl());
			strata.println("role_partial_stratum\tcopies\tunique_sequences\tcovered_copies\tcovered_unique\tuncovered_copies\tuncovered_unique");
			for(String key:sorted(totals)){final Composition total=totals.get(key),yes=coveredGroups.get(key);final long cp=yes==null?0:yes.copies;final int uq=yes==null?0:yes.unique.size();strata.print(b.clear().append(key).tab().append(total.copies).tab().append(total.unique.size()).tab().append(cp).tab().append(uq).tab().append(total.copies-cp).tab().append(total.unique.size()-uq).nl());}
		}finally{RrnaResourceIO.finish(steps,fasta,rows,composition,summary,strata);}
		System.err.println("MITO_COVERAGE_SELECTION_PASS copies="+members.size()+" unique="+unique.size()+" selected="+selected.length+" uncovered="+(unique.size()-covered)+" cut="+cut+" mode="+mode+" canonical_models="+canonicalModels+" actual_alignment_calls="+alignmentCalls+" matrix_rows="+((long)measuredUnique*models.size()));
	}
	/** Frozen model-order ties; no cap. Rows are primitive bitsets of covered unique sequences. */
	public static int[] greedy(long[][] cover,int count){
		require(count>0 && cover.length>0,"Nonempty coverage universe and candidate pool required");final int words=(count+63)/64;
		final long[] remaining=new long[words];Arrays.fill(remaining,-1L);if((count&63)!=0){remaining[words-1]=(1L<<(count&63))-1;}
		for(long[] bits:cover){require(bits.length==words,"Coverage bitset length must match universe");}
		final IntList chosen=new IntList();while(true){int best=-1,gain=0;for(int m=0;m<cover.length;m++){int n=0;for(int w=0;w<words;w++){n+=Long.bitCount(cover[m][w]&remaining[w]);}if(n>gain){gain=n;best=m;}}
			if(best<0){break;}chosen.add(best);for(int w=0;w<words;w++){remaining[w]&=~cover[best][w];}}
		return chosen.toArray();
	}
	private static void writeComposition(ByteStreamWriter out,ByteBuilder b,Model model,String scope,ObjectMap<String,Composition> groups){
		for(String name:sorted(groups)){final Composition c=groups.get(name);out.print(b.clear().append(model.key).tab().append(scope).tab().append(name).tab().append(c.copies).tab().append(c.unique.size()).nl());}
	}
	private static void composition(ObjectMap<String,Composition> groups,String phylum,String uid){Composition c=groups.get(phylum);if(c==null){c=new Composition();groups.put(phylum,c);}c.copies++;c.unique.add(uid);}
	private String phylum(int taxid){int i=lineageIds.get(taxid);if(i<0){i=lineages.size();lineages.add(NcrnaMitoClusterEvaluation.lineage(tree,taxid,TaxTree.PHYLUM_E));lineageIds.put(taxid,i);}return lineages.get(i);}
	private static <V> ArrayList<String> sorted(ObjectMap<String,V> map){final ArrayList<String> a=new ArrayList<String>();for(String s:map.keys()){if(s!=null){a.add(s);}}Collections.sort(a);return a;}
	private static void validate(Sequence s,Model m){require(Float.isFinite(s.identity[m.index]) && s.identity[m.index]>=0 && s.identity[m.index]<=1 && Float.isFinite(s.trace[m.index]) && s.trace[m.index]>=0 && s.trace[m.index]<=1 && s.paired[m.index]>=0 && s.paired[m.index]<=Math.min(s.bases.length,m.bases.length),"Invalid aligned-pair evidence for "+s.id+"/"+m.key);}
	@Override public void accumulate(Worker w){ok&=w.ok;alignmentCalls+=w.alignmentCalls;if(w.failure!=null){w.failure.printStackTrace(System.err);}}
	@Override public boolean success(){return ok;}
	@Override public ReadWriteLock rwlock(){return null;}
	static final class Worker extends Thread{
		Worker(NcrnaMitoCoverageSelect owner_){owner=owner_;setDaemon(true);}
		@Override public void run(){try{final AlignmentStats stats=new AlignmentStats(true);
			for(int i;!owner.abort && (i=owner.next.getAndIncrement())<owner.unique.size();){final Sequence s=owner.unique.get(i);if(!s.measured){continue;}
				for(Model m:owner.models){if(m.canonical!=m.index){s.identity[m.index]=s.identity[m.canonical];s.trace[m.index]=s.trace[m.canonical];s.paired[m.index]=s.paired[m.canonical];continue;}
					stats.clear();stats.doTrace=true;final boolean shorter=s.bases.length<=m.bases.length;final byte[] q=shorter?s.bases:m.bases,r=shorter?m.bases:s.bases;
					alignmentCalls++;s.identity[m.index]=QuantumAligner.alignAndTraceStatic(q,r,stats);s.trace[m.index]=stats.identity;s.paired[m.index]=stats.matches+stats.subs+stats.ns;
					require(s.paired[m.index]+stats.ins==q.length && s.paired[m.index]+stats.dels==stats.rStop-stats.rStart+1 && stats.rStart>=0 && stats.rStop<r.length,"Native shorter-query trace must conserve both sequences");validate(s,m);
				}
			}ok=!owner.abort;
		}catch(Throwable t){failure=t;owner.abort=true;}}
		final NcrnaMitoCoverageSelect owner;boolean ok;Throwable failure;long alignmentCalls;
	}
	private static final class Key{Key(byte[] b){bases=b;hash=Arrays.hashCode(b);}final byte[] bases;final int hash;@Override public int hashCode(){return hash;}@Override public boolean equals(Object x){return x instanceof Key && Arrays.equals(bases,((Key)x).bases);}}
	private static final class Sequence{Sequence(byte[] b,int i,int models){bases=b;index=i;id="U"+i;identity=new float[models];trace=new float[models];paired=new int[models];}final byte[] bases;final int index;final String id;final float[] identity,trace;final int[] paired;boolean construction,panel,measured;}
	private static final class Model{Model(String k,String r,String t,String h,byte[] b,int i){key=k;role=r;header=h;bases=b;index=i;output="mitochondria_coverage_"+r+'_'+t;}final String key,role,header,output;final byte[] bases;final int index;int canonical;final ObjectMap<String,Composition> sourceComposition=new ObjectMap<String,Composition>(String.class,Composition.class);}
	private static final class Source{Source(String i,String u,int t,int l,boolean p,boolean ps,String r,String s){id=i;uid=u;taxid=t;length=l;partial=p;pseudo=ps;role=r;path=s;require(t>0 && l>0,"Source member requires taxid and length");}final String id,uid,role,path;final int taxid,length;final boolean partial,pseudo;String model="NA";}
	private static final class Member{Member(String r,String i,Sequence s,int t,boolean p,boolean ps,String m,String str,String source_){role=r;id=i;sequence=s;taxid=t;partial=p;pseudo=ps;sourceModel=m;stratum=str;source=source_;}final String role,id,sourceModel,stratum,source;final Sequence sequence;final int taxid;final boolean partial,pseudo;}
	private static final class Panel{Panel(LineParser1 p){id=p.parseString(0);assembly=p.parseString(2);contig=p.parseString(3);start=p.parseInt(4);stop=p.parseInt(5);strand=p.parseString(6);taxid=p.parseInt(37);final String a=p.parseString(29);require(a.startsWith("LOCUS:") && start>0 && stop>=start && taxid>0,"Exact panel provenance required");truthId=a.substring(6);}final String id,assembly,contig,strand,truthId;final int start,stop,taxid;boolean bound,seen,partial;String geometry;}
	private static final class Composition{long copies;final ObjectSet<String> unique=new ObjectSet<String>(String.class);}
	private static boolean bool(String s){require(s.equals("true") || s.equals("false"),"Literal source boolean required");return s.equals("true");}
	private static boolean bit(LineParser1 p,int n){require(p.termEquals("0",n) || p.termEquals("1",n),"Frozen truth bit must be explicit");return p.termEquals("1",n);}
	private String required(String k){final String v=flags.get(k);require(v!=null,"Missing flag "+k);return v;}
	private String value(String k,String fallback){final String v=flags.get(k);return v==null?fallback:v;}
	private int integer(String k,int fallback){return Parse.parseIntKMG(value(k,Integer.toString(fallback)));}
	private static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final String[] INPUTS={"corpus","ssumodels","nohitmodels","ssumembers","nohitmembers","panelcases","panelgenes","truth","tree"};
	private static final String KNOWN="|corpus|ssumodels|nohitmodels|ssumembers|nohitmembers|panelcases|panelgenes|truth|tree|verified|outprefix|mode|matrix|expected|expectedpanel|expectedssu|expectednohit|expectedssucopies|expectednohitcopies|pilotsize|t|cut|";
	private final ObjectMap<String,String> flags=new ObjectMap<String,String>(String.class,String.class);private ObjectMap<String,String> pins;private TaxTree tree;
	private String prefix,mode;private int corpusExpected,panelExpected,ssuExpected,nohitExpected,pilotSize,threads,canonicalModels,measuredUnique;private float cut;private long alignmentCalls;
	private final ArrayList<Model> models=new ArrayList<Model>();private final ObjectMap<String,Model> modelByKey=new ObjectMap<String,Model>(String.class,Model.class);
	private final ObjectMap<String,Source> sourceMembers=new ObjectMap<String,Source>(String.class,Source.class);
	private final ObjectMap<String,Sequence> sourceSequences=new ObjectMap<String,Sequence>(String.class,Sequence.class);
	private final IntHashMap lineageIds=new IntHashMap();private final ArrayList<String> lineages=new ArrayList<String>();
	private final ArrayList<Member> members=new ArrayList<Member>();private final ArrayList<Sequence> unique=new ArrayList<Sequence>();private final ObjectMap<Key,Sequence> cache=new ObjectMap<Key,Sequence>(Key.class,Sequence.class);
	private final AtomicInteger next=new AtomicInteger();private volatile boolean abort;private boolean ok=true;
}
