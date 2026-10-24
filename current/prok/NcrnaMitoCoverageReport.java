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

/** Reports saved-matrix selections without aligning or changing the selector.
 * Copies and globally unique sequences are conserved independently. Per-lineage
 * unique counts are incidences: their sum need not equal a global unique count.
 * @author Brian, Ganyu */
public final class NcrnaMitoCoverageReport {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args,NcrnaMitoCoverageReport.class,false);
		try{final NcrnaMitoCoverageReport x=new NcrnaMitoCoverageReport();x.parse(pp.args);x.run();}
		finally{Shared.closeStream(pp.outstream);}
	}
	private void parse(String[] args){
		for(String arg:args){final String[] p=Parse.splitOnFirst(arg,'=');require(p.length==2,"Use flag=value");final String key=p[0].toLowerCase(Locale.ROOT);
			require(KNOWN.contains("|"+key+"|") && flags.put(key,p[1])==null,"Unknown/repeated flag: "+key);}
		for(String key:new String[]{"runs","verified","outprefix"}){require(flags.get(key)!=null,"Missing report binding: "+key);}
		final String[] values=value("cuts",".60,.62,.64,.66,.68,.70").split(",");cuts=new float[values.length];boolean found=false;target=Float.parseFloat(value("target",".66"));
		for(int i=0;i<cuts.length;i++){cuts[i]=Float.parseFloat(values[i]);require(Float.isFinite(cuts[i]) && cuts[i]>0 && cuts[i]<=1 && (i==0 || cuts[i]>cuts[i-1]),"Cuts must be unique increasing fractions");found|=cuts[i]==target;}require(found,"Target must be one of the declared cuts");
		pins=EukPanelGffValidation.hashes(flags.get("verified"));bound(flags.get("runs"));Shared.setThreads(1);ReadWrite.setZipThreads(1);
	}
	private void run(){
		final String out=flags.get("outprefix");RrnaResourceIO.fresh(out,".curve.tsv",".residual_lineages.tsv",".models.tsv",".inputs.tsv");
		final ArrayList<Run> runs=new ArrayList<Run>();read(flags.get("runs"),"cut\tprefix",2,p->{final int i=runs.size();require(i<cuts.length && p.parseFloat(0)==cuts[i],"Roster must contain each declared cut once in increasing order");final Run r=new Run(cuts[i],p.parseString(1));for(String suffix:SUFFIXES){bound(r.prefix+suffix);}runs.add(r);});
		require(runs.size()==cuts.length,"Every requested threshold must be present");
		for(Run r:runs){loadSummary(r);loadModels(r);loadSteps(r);loadMembers(r);checkStrata(r);loadComposition(r);}
		write(runs,out);
	}
	private void loadSummary(Run r){
		read(r.prefix+".summary.tsv",SUMMARY_HEADER,16,p->{require(r.copies<0,"Exactly one summary data row required");r.copies=p.parseLong(0);r.unique=p.parseInt(1);r.modelCount=p.parseInt(2);r.selected=p.parseInt(3);r.coveredUnique=p.parseInt(5);r.coveredCopies=p.parseLong(7);r.partial=p.parseLong(9);r.panel=p.parseLong(10);r.canonical=p.parseInt(13);r.matrixRows=p.parseLong(15);
			require(r.copies>0 && r.unique>0 && r.unique<=r.copies && r.modelCount>0 && r.selected>=0 && r.selected<=r.modelCount && r.coveredUnique>=0 && r.coveredCopies>=0 && r.coveredUnique+p.parseInt(6)==r.unique && r.coveredCopies+p.parseLong(8)==r.copies,"Summary copy and unique denominators must conserve all members");
			require(p.parseFloat(4)==r.cut && p.termEquals("SHORTER_QUERY_COVERAGE_NOT_FULL_MODEL_CALLER_VERIFICATION",11) && p.termEquals("select",12) && p.parseLong(14)==0 && r.canonical>0 && r.canonical<=r.modelCount && r.matrixRows==(long)r.unique*r.modelCount,"Report requires zero-alignment saved full-matrix selection with explicit shorter-query semantics");});require(r.copies>0,"Missing summary row");
	}
	private void loadModels(Run r){
		final String hash=bound(r.prefix+".models.tsv");if(modelHash==null){modelHash=hash;}require(hash.equals(modelHash),"All thresholds must use the identical model roster");
		read(r.prefix+".models.tsv",MODEL_HEADER,6,p->{final Model m=new Model(p);require(r.models.put(m.id,m)==null,"Unique model IDs required");});
		int canonical=0;for(Model m:r.models.values()){if(m!=null){final Model c=r.models.get(m.canonical);require(c!=null && c.id.equals(c.canonical),"Every alias must bind a canonical model");if(m.id.equals(m.canonical)){canonical++;}}}
		require(r.models.size()==r.modelCount && canonical==r.canonical,"Model and canonical counts must match the full matrix");
	}
	private void loadSteps(Run r){
		final long[] total={0,0};read(r.prefix+".selection.tsv",STEP_HEADER,6,p->{final int step=p.parseInt(0);final Model m=r.models.get(p.parseString(1));
			require(step==r.steps.size()+1 && m!=null && m.step==0 && m.canonical.equals(m.id),"Selection order must contain unique canonical models");m.step=step;m.gain=p.parseInt(2);m.copies=p.parseLong(3);require(m.gain>0 && m.copies>=m.gain,"Each greedy step must add positive unique coverage");total[0]+=m.gain;total[1]+=m.copies;
			require(p.parseLong(4)==total[0] && p.parseLong(5)==r.unique-total[0],"Selection cumulative and residual counts must conserve the universe");r.steps.add(m);});
		require(r.steps.size()==r.selected && total[0]==r.coveredUnique && total[1]==r.coveredCopies,"Greedy step totals must agree with the summary");
	}
	private void loadMembers(Run r){
		final boolean first=baseline.isEmpty();final ObjectSet<String> seen=new ObjectSet<String>(String.class);final ObjectMap<String,Unique> uniques=new ObjectMap<String,Unique>(String.class,Unique.class);final long[] counters={0,0};
		read(r.prefix+".members.tsv.gz",MEMBER_HEADER,19,p->{final String role=p.parseString(0),id=p.parseString(1),uid=p.parseString(2),phylum=p.parseString(4),partial=p.parseString(5),assigned=p.parseString(9),best=p.parseString(10),pool=p.parseString(15);
			require((role.equals("CONSTRUCTION") || role.equals("PANEL_INCLUDED")) && (partial.equals("true") || partial.equals("false")) && p.parseInt(3)>0,"Explicit role, partial state and taxonomy required");
			final String key=role+'\t'+id;require(seen.add(key),"Each original copy must occur once per threshold");
			final String signature=uid+'\t'+p.parseString(3)+'\t'+phylum+'\t'+partial+'\t'+p.parseString(6)+'\t'+p.parseString(7)+'\t'+p.parseString(8)+'\t'+p.parseString(18);
			final float poolIdentity=p.parseFloat(16);require(Float.isFinite(poolIdentity) && poolIdentity>=0 && poolIdentity<=1 && r.models.get(pool)!=null,"Best-pool evidence must name a real model and finite identity");
			if(first){baseline.put(key,signature);baselinePool.put(uid,pool+'\t'+Float.toString(poolIdentity));}
			else{require(signature.equals(baseline.get(key)) && (pool+'\t'+Float.toString(poolIdentity)).equals(baselinePool.get(uid)),"Thresholds must retain identical member provenance and underlying pool evidence");}
			final boolean covered=!assigned.equals("NA");final Model selected=r.models.get(assigned);final Model winner=r.models.get(best);float identity=-1;
			if(!best.equals("NA")){identity=p.parseFloat(11);require(winner!=null && winner.step>0 && Float.isFinite(identity) && identity>=0 && identity<=1 && identity<=poolIdentity,"Selected winner must belong to selected models and not exceed pool best");}
			else{require(r.selected==0 && p.termEquals("NA",11),"Missing selected winner requires an empty selection");}
			require(covered==(identity>=r.cut) && (covered?selected!=null && selected.step>0:p.termEquals("NA",9)) && p.termEquals(covered?"COVERED_SHORTER_QUERY":"UNCOVERED_AT_CUT",17),"Assigned coverage must agree with selected identity and declared threshold");
			final Unique u=new Unique(assigned,best,identity,pool,poolIdentity);final Unique prior=uniques.get(uid);if(prior==null){uniques.put(uid,u);}else{require(prior.same(u),"Exact-sequence copies cannot have different selection evidence");}
			r.all.add(uid,covered);group(r.strata,role+":"+(partial.equals("true")?"PARTIAL":"NOT_MARKED_PARTIAL")).add(uid,covered);
			group(r.lineages,phylum).add(uid,covered);group(r.detail,role+'\t'+partial+'\t'+phylum).add(uid,covered);
			if(covered){selected.assigned.add(uid,true);group(selected.newGroups,phylum).add(uid,true);}
			if(partial.equals("true")){counters[0]++;}if(role.equals("PANEL_INCLUDED")){counters[1]++;}
		});
		require(seen.size()==baseline.size() && r.all.copies==r.copies && r.all.unique.size()==r.unique && r.all.covered==r.coveredCopies && r.all.yes.size()==r.coveredUnique && counters[0]==r.partial && counters[1]==r.panel,"Member-derived coverage must independently reproduce all summary denominators");
		for(Model m:r.steps){require(m.assigned.copies==m.copies && m.assigned.unique.size()==m.gain,"Member assignments must reproduce each greedy step gain: "+m.id);}
	}
	private void checkStrata(Run r){
		final ObjectSet<String> seen=new ObjectSet<String>(String.class);read(r.prefix+".strata.tsv",STRATA_HEADER,7,p->{final String key=p.parseString(0);final Group g=r.strata.get(key);require(g!=null && seen.add(key),"Every reported stratum must uniquely match member rows");
			require(g.copies==p.parseLong(1) && g.unique.size()==p.parseInt(2) && g.covered==p.parseLong(3) && g.yes.size()==p.parseInt(4) && g.copies-g.covered==p.parseLong(5) && g.unique.size()-g.yes.size()==p.parseInt(6),"Stratum counts must be derived from the same full membership");});require(seen.size()==r.strata.size(),"Missing member-derived stratum");
	}
	private void loadComposition(Run r){
		final ObjectSet<String> seen=new ObjectSet<String>(String.class);read(r.prefix+".composition.tsv",COMPOSITION_HEADER,5,p->{final String model=p.parseString(0),scope=p.parseString(1),phylum=p.parseString(2);final Model m=r.models.get(model);require(m!=null && seen.add(model+'\t'+scope+'\t'+phylum),"Composition rows require unique model/scope/lineage keys");
			final Count n=new Count(p.parseLong(3),p.parseInt(4));require(n.copies>0 && n.unique>0 && n.unique<=n.copies,"Per-lineage unique counts cannot exceed copy counts");
			if(scope.equals("NEWLY_COVERED_TARGETS")){final Group g=m.newGroups.get(phylum);require(m.step>0 && g!=null && g.copies==n.copies && g.unique.size()==n.unique,"Newly-covered composition must match the actual assigned member rows");m.newSeen++;}
			else{require(scope.equals("SOURCE_CLUSTER_FOR_SELECTED_"+m.canonical) && r.models.get(m.canonical).step>0,"Source composition must explicitly bind a selected canonical model or its alias");require(m.source.put(phylum,n)==null,"Duplicate original-cluster lineage");}
		});for(Model m:r.steps){require(m.newSeen==m.newGroups.size(),"Missing newly-covered lineage row for selected model: "+m.id);}
	}
	private void write(ArrayList<Run> runs,String prefix){
		final ByteStreamWriter curve=RrnaResourceIO.writer(prefix+".curve.tsv"),residual=RrnaResourceIO.writer(prefix+".residual_lineages.tsv"),models=RrnaResourceIO.writer(prefix+".models.tsv"),inputs=RrnaResourceIO.writer(prefix+".inputs.tsv");final ByteBuilder b=new ByteBuilder();
		try{curve.println("cut\tstratum\tselected_models\tcopies\tunique_sequences\tcovered_copies\tcovered_unique\tuncovered_copies\tuncovered_unique\tcopy_coverage\tunique_coverage\tnew_alignment_calls\tmatrix_rows\tinterpretation");
			residual.println("cut\tscope\trole\tpartial\tphylum\tcopies\tunique_sequences\tcovered_copies\tcovered_unique\tuncovered_copies\tuncovered_unique\tunique_counts_additive_across_lineages");
			models.println("cut\tstep\tselected_model\tevidence_model\toutput_id\tmodel_length\tscope\tcopies\tglobal_unique_sequences\tlineage_unique_incidence_sum\tdominant_lineages_by_copies\tcopy_tie_count\tmax_lineage_copies\tdominant_lineages_by_unique\tunique_tie_count\tmax_lineage_unique\toriginal_header");inputs.println("cut\trole\tpath\tsha80");
			inputs.print(b.clear().append("NA\truns\t").append(flags.get("runs")).tab().append(bound(flags.get("runs"))).nl());
			for(Run r:runs){curve(curve,b,r,"ALL",r.all);for(String key:keys(r.strata)){curve(curve,b,r,key,r.strata.get(key));}
				for(String suffix:SUFFIXES){inputs.print(b.clear().appendSlow(r.cut).tab().append(suffix).tab().append(r.prefix+suffix).tab().append(bound(r.prefix+suffix)).nl());}
				if(r.cut!=target){continue;}residual(residual,b,r,"ALL","ALL\tALL\tALL",r.all);
				for(String phylum:keys(r.lineages)){residual(residual,b,r,"PHYLUM","ALL\tALL\t"+phylum,r.lineages.get(phylum));}
				for(String key:keys(r.detail)){residual(residual,b,r,"ROLE_PARTIAL_PHYLUM",key,r.detail.get(key));}
				for(Model m:r.steps){for(String id:keys(r.models)){final Model alias=r.models.get(id);if(alias.canonical.equals(m.id)){model(models,b,r,m,alias,"SOURCE_CLUSTER",alias.source,-1);}}
					final ObjectMap<String,Count> counts=new ObjectMap<String,Count>(String.class,Count.class);for(String lineage:keys(m.newGroups)){final Group g=m.newGroups.get(lineage);counts.put(lineage,new Count(g.copies,g.unique.size()));}model(models,b,r,m,m,"NEWLY_COVERED_TARGETS",counts,m.gain);}
			}
		}finally{RrnaResourceIO.finish(curve,residual,models,inputs);}
		System.err.println("MITO_COVERAGE_REPORT_PASS cuts="+runs.size()+" copies="+baseline.size()+" target="+target+" new_alignments=0");
	}
	private static void curve(ByteStreamWriter out,ByteBuilder b,Run r,String scope,Group g){out.print(b.clear().appendSlow(r.cut).tab().append(scope).tab().append(r.selected).tab().append(g.copies).tab().append(g.unique.size()).tab().append(g.covered).tab().append(g.yes.size()).tab().append(g.copies-g.covered).tab().append(g.unique.size()-g.yes.size()).tab().append(g.covered/(double)g.copies,9).tab().append(g.yes.size()/(double)g.unique.size(),9).append("\t0\t").append(r.matrixRows).append("\tSHORTER_QUERY_COVERAGE_NOT_CALLER_RECALL\n"));}
	private static void residual(ByteStreamWriter out,ByteBuilder b,Run r,String scope,String key,Group g){out.print(b.clear().appendSlow(r.cut).tab().append(scope).tab().append(key).tab().append(g.copies).tab().append(g.unique.size()).tab().append(g.covered).tab().append(g.yes.size()).tab().append(g.copies-g.covered).tab().append(g.unique.size()-g.yes.size()).append("\tfalse\n"));}
	private static void model(ByteStreamWriter out,ByteBuilder b,Run r,Model selected,Model evidence,String scope,ObjectMap<String,Count> counts,int globalUnique){
		long copies=0,incidences=0,maxCopies=0;int maxUnique=0,copyTies=0,uniqueTies=0;final StringBuilder copyNames=new StringBuilder(),uniqueNames=new StringBuilder();
		for(String name:keys(counts)){final Count c=counts.get(name);copies+=c.copies;incidences+=c.unique;maxCopies=Math.max(maxCopies,c.copies);maxUnique=Math.max(maxUnique,c.unique);}
		for(String name:keys(counts)){final Count c=counts.get(name);if(c.copies==maxCopies){if(copyTies++>0){copyNames.append(';');}copyNames.append(name);}if(c.unique==maxUnique){if(uniqueTies++>0){uniqueNames.append(';');}uniqueNames.append(name);}}
		b.clear().appendSlow(r.cut).tab().append(selected.step).tab().append(selected.id).tab().append(evidence.id).tab().append(selected.output).tab().append(evidence.length).tab().append(scope).tab().append(copies).tab();if(globalUnique<0){b.append("NA");}else{b.append(globalUnique);}
		out.print(b.tab().append(incidences).tab().append(copyTies==0?"NA":copyNames.toString()).tab().append(copyTies).tab().append(maxCopies).tab().append(uniqueTies==0?"NA":uniqueNames.toString()).tab().append(uniqueTies).tab().append(maxUnique).tab().append(evidence.header).nl());
	}
	private String bound(String path){final String hash=pins.get(path);require(hash!=null,"Missing producer-verified input: "+path);return hash;}
	private String value(String key,String fallback){final String s=flags.get(key);return s==null?fallback:s;}
	private interface Row {void accept(LineParser1 p);}
	private static void read(String path,String header,int terms,Row row){final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');try{RrnaResourceIO.header(in,header,path);for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==terms,"Unexpected schema: "+path);row.accept(p);}}finally{require(!in.close(),"Table read failed: "+path);}}
	private static Group group(ObjectMap<String,Group> map,String key){Group g=map.get(key);if(g==null){g=new Group();map.put(key,g);}return g;}
	private static <V> ArrayList<String> keys(ObjectMap<String,V> map){final ArrayList<String> a=new ArrayList<String>();for(String key:map.keys()){if(key!=null){a.add(key);}}Collections.sort(a);return a;}
	private static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final class Group{long copies,covered;final ObjectSet<String> unique=new ObjectSet<String>(String.class),yes=new ObjectSet<String>(String.class);void add(String uid,boolean hit){copies++;unique.add(uid);if(hit){covered++;yes.add(uid);}}}
	private static final class Count{Count(long c,int u){copies=c;unique=u;}final long copies;final int unique;}
	private static final class Unique{Unique(String a,String b,float i,String p,float pi){assigned=a;best=b;identity=i;pool=p;poolIdentity=pi;}boolean same(Unique u){return assigned.equals(u.assigned) && best.equals(u.best) && identity==u.identity && pool.equals(u.pool) && poolIdentity==u.poolIdentity;}final String assigned,best,pool;final float identity,poolIdentity;}
	private static final class Model{Model(LineParser1 p){id=p.parseString(0);output=p.parseString(1);length=p.parseInt(2);canonical=p.parseString(3);header=p.parseString(5);require(length>0,"Selected models need positive sequence lengths");}final String id,output,canonical,header;final int length;int step,gain,newSeen;long copies;final Group assigned=new Group();final ObjectMap<String,Group> newGroups=new ObjectMap<String,Group>(String.class,Group.class);final ObjectMap<String,Count> source=new ObjectMap<String,Count>(String.class,Count.class);}
	private static final class Run{Run(float c,String p){cut=c;prefix=p;}final float cut;final String prefix;long copies=-1,coveredCopies,partial,panel,matrixRows;int unique,modelCount,canonical,selected,coveredUnique;final Group all=new Group();final ObjectMap<String,Model> models=new ObjectMap<String,Model>(String.class,Model.class);final ArrayList<Model> steps=new ArrayList<Model>();final ObjectMap<String,Group> strata=new ObjectMap<String,Group>(String.class,Group.class),lineages=new ObjectMap<String,Group>(String.class,Group.class),detail=new ObjectMap<String,Group>(String.class,Group.class);}
	private static final String KNOWN="|runs|verified|outprefix|cuts|target|";
	private static final String[] SUFFIXES={".members.tsv.gz",".models.tsv",".composition.tsv",".summary.tsv",".strata.tsv",".selection.tsv"};
	public static final String SUMMARY_HEADER="copies\tunique_sequences\tmodels\tselected_models\tcut\tcovered_unique\tuncovered_unique\tcovered_copies\tuncovered_copies\tpartial_copies\tpanel_included_copies\tinterpretation\tmode\tcanonical_models\tactual_alignment_calls\tmatrix_rows";
	public static final String MODEL_HEADER="model_id\toutput_id\tlength\tcanonical_model\tsource_namespace\toriginal_header";
	public static final String STEP_HEADER="step\tmodel_id\tnew_unique\tnew_copies\tcovered_unique\tremaining_unique";
	public static final String MEMBER_HEADER="role\tlocus_id\tunique_id\ttaxid\tphylum\tpartial\tpseudo\tstratum\tsource_model\tassigned_cover_model\tbest_selected_model\tselected_identity\tselected_trace_identity\tmember_paired_fraction\tmodel_paired_fraction\tbest_pool_model\tbest_pool_identity\tstatus\tprovenance";
	public static final String STRATA_HEADER="role_partial_stratum\tcopies\tunique_sequences\tcovered_copies\tcovered_unique\tuncovered_copies\tuncovered_unique";
	public static final String COMPOSITION_HEADER="model_id\tscope\tphylum\tcopies\tunique_sequences";
	private final ObjectMap<String,String> flags=new ObjectMap<String,String>(String.class,String.class),baseline=new ObjectMap<String,String>(String.class,String.class),baselinePool=new ObjectMap<String,String>(String.class,String.class);private ObjectMap<String,String> pins;private float[] cuts;private float target;private String modelHash;
}
