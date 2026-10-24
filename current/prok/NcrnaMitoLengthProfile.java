package prok;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import map.ObjectMap;
import map.ObjectIntMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;
import structures.IntList;

/** Derives per-model length bounds from every covered member of a saved matrix.
 * A member contributes to every selected model meeting the cut, irrespective of
 * its greedy assignment. No alignment or sequence/model selection is performed.
 * @author Brian, Ganyu */
public final class NcrnaMitoLengthProfile {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args,NcrnaMitoLengthProfile.class,false);
		try{final NcrnaMitoLengthProfile x=new NcrnaMitoLengthProfile();x.parse(pp.args);x.run();}
		finally{Shared.closeStream(pp.outstream);}
	}
	private void parse(String[] args){
		for(String arg:args){final String[] p=Parse.splitOnFirst(arg,'=');require(p.length==2,"Use flag=value");final String key=p[0].toLowerCase(Locale.ROOT);
			require(KNOWN.contains("|"+key+"|") && flags.put(key,p[1])==null,"Unknown/repeated flag: "+key);}
		cut=Float.parseFloat(value("cut",".66"));expectedSelected=Parse.parseIntKMG(value("expectedselected","153"));
		if(flags.contains("maxwindow")){maxWindow=Parse.parseIntKMG(required("maxwindow"));require(maxWindow>0,"Explicit maxwindow must be positive");}
		require(Float.isFinite(cut) && cut>0 && cut<=1 && expectedSelected>0,"A finite coverage fraction and positive selected-model count are required");
		pins=EukPanelGffValidation.hashes(required("verified"));for(String key:BINDINGS){bound(key);}for(String key:INPUTS){bound(key);}
		prefix=required("outprefix");RrnaResourceIO.fresh(prefix,".bounds.tsv",".profile.tsv",".summary.tsv",".inputs.tsv");
		if(maxWindow>0){RrnaResourceIO.fresh(prefix,".window_summary.tsv",".window_lineages.tsv");}
		Shared.setThreads(1);ReadWrite.setZipThreads(1);
	}
	private void run(){
		loadModels();loadSelection();loadUnique();loadMembers();loadSummary();checkSelected();loadMatrix();write();
	}
	private void loadModels(){
		final ObjectSet<String> outputs=new ObjectSet<String>(String.class);
		read(required("models"),NcrnaMitoCoverageReport.MODEL_HEADER,6,p->{final Model m=new Model(p);
			require(byModel.put(m.id,m)==null && outputs.add(m.output),"Model and output IDs must be unique");
			if(m.id.equals(m.canonical)){m.original=m;canonicalModels++;}else{m.original=byModel.get(m.canonical);require(m.original!=null && m.original.original==m.original && m.original.length==m.length,"Aliases must follow their exact canonical model with the same length");}
			models.add(m);});require(!models.isEmpty(),"A full matrix model roster is required");
	}
	private void loadSelection(){
		final long[] cumulative={0};read(required("selection"),NcrnaMitoCoverageReport.STEP_HEADER,6,p->{final Model m=byModel.get(p.parseString(1));
			require(p.parseInt(0)==selected.size()+1 && m!=null && m.original==m && m.step==0,"Selection must contain unique canonical models in step order");
			m.step=p.parseInt(0);m.newUnique=p.parseInt(2);m.newCopies=p.parseLong(3);require(m.newUnique>0 && m.newCopies>=m.newUnique,"Every selected model must have a positive saved greedy gain");
			cumulative[0]+=m.newUnique;require(p.parseLong(4)==cumulative[0] && p.parseLong(5)>=0,"Saved greedy cumulative coverage must conserve unique members");
			final long universe=p.parseLong(4)+p.parseLong(5);if(selectionUniverse<0){selectionUniverse=universe;}require(universe==selectionUniverse,"Selection steps must use one fixed universe");selected.add(m);
		});require(selected.size()==expectedSelected,"Selected-model count differs from the explicitly requested library");
		if(maxWindow>0){for(Model m:selected){m.identityMinimum=MitoSsuModelPolicy.minimumLength(m.length,cut);if(m.identityMinimum>maxWindow){m.windowCovered=new BitSet();ineligibleModels++;}}}
	}
	private void loadUnique(){
		RrnaResourceIO.stream(required("unique"),r->{final String id=RrnaResourceIO.first(r.id);require(id.equals("U"+uniques.size()) && r.length()>0,"Unique FASTA must retain complete ordered matrix IDs and positive lengths");
			final Unique u=new Unique(id,r.length(),uniques.size());uniques.add(u);byUnique.put(id,u);});
		require(!uniques.isEmpty() && uniques.size()==selectionUniverse,"Unique FASTA and saved selection must bind the same full universe");
	}
	private void loadMembers(){
		final ObjectSet<String> seen=new ObjectSet<String>(String.class);
		read(required("members"),NcrnaMitoCoverageReport.MEMBER_HEADER,19,p->{
			final String role=p.parseString(0);require(role.equals("CONSTRUCTION") || role.equals("PANEL_INCLUDED"),"Keep construction and included-panel members explicit");
			require(seen.add(role+'\t'+p.parseString(1)),"An original member copy must not be counted twice");
			final Unique u=byUnique.get(p.parseString(2));require(u!=null,"Every member must bind a sequence in the exact unique FASTA");
			final boolean partial=bit(p,5);final String assigned=p.parseString(9),best=p.parseString(10);final Model a=byModel.get(assigned),b=byModel.get(best);final float identity=p.parseFloat(11);
			require(b!=null && b.step>0 && Float.isFinite(identity) && identity>=0 && identity<=1,"Best selected evidence must name a selected model and finite identity");
			require((identity>=cut)==!assigned.equals("NA") && (assigned.equals("NA") || a!=null && a.step>0) && p.termEquals(identity>=cut?"COVERED_SHORTER_QUERY":"UNCOVERED_AT_CUT",17),"Member assignment/status must agree with the requested coverage cut");
			if(u.copies==0){u.assigned=a;u.best=b;u.bestIdentity=identity;}else{require(u.assigned==a && u.best==b && u.bestIdentity==identity,"Exact-sequence copies must carry identical selected-model evidence");}
			u.copies++;copies++;if(partial){partialCopies++;}if(role.equals("PANEL_INCLUDED")){panelCopies++;}
			if(maxWindow>0){final String lineage=p.parseString(4);require(!lineage.isEmpty(),"Window analysis requires explicit per-copy lineage metadata");
				if(!lineageIndex.contains(lineage)){lineageIndex.put(lineage,lineageNames.size());lineageNames.add(lineage);}
				copyUnique.add(u.index);copyLineage.add(lineageIndex.get(lineage));}
			if(a!=null){a.assignedCopies++;if(u.copies==1){a.assignedUnique++;}}
		});
		for(Unique u:uniques){require(u.copies>0,"Unique FASTA may not contain members absent from the full copy roster: "+u.id);}
		for(Model m:selected){require(m.assignedUnique==m.newUnique && m.assignedCopies==m.newCopies,"Saved greedy gains must match copy/unique member assignments: "+m.id);}
	}
	private void loadSummary(){
		final int[] rows={0};read(required("summary"),NcrnaMitoCoverageReport.SUMMARY_HEADER,16,p->{require(++rows[0]==1,"Exactly one saved selection summary is required");
			long coveredCopies=0;int coveredUnique=0;for(Model m:selected){coveredCopies+=m.newCopies;coveredUnique+=m.newUnique;}
			require(p.parseLong(0)==copies && p.parseInt(1)==uniques.size() && p.parseInt(2)==models.size() && p.parseInt(3)==selected.size() && p.parseFloat(4)==cut
				&& p.parseInt(5)==coveredUnique && p.parseInt(6)==uniques.size()-coveredUnique && p.parseLong(7)==coveredCopies && p.parseLong(8)==copies-coveredCopies
				&& p.parseLong(9)==partialCopies && p.parseLong(10)==panelCopies,"Summary must exactly conserve bound members, selected models and strata");
			require(p.termEquals("SHORTER_QUERY_COVERAGE_NOT_FULL_MODEL_CALLER_VERIFICATION",11) && p.termEquals("select",12) && p.parseInt(13)==canonicalModels && p.parseLong(14)==0 && p.parseLong(15)==(long)uniques.size()*models.size(),"Require saved full-matrix selection, not a pilot or new alignment run");
		});require(rows[0]==1,"Missing saved selection summary");
	}
	private void checkSelected(){
		final int[] n={0};RrnaResourceIO.stream(required("selected"),r->{require(n[0]<selected.size(),"Extra selected FASTA model");final Model m=selected.get(n[0]++);
			require(RrnaResourceIO.first(r.id).equals(m.output) && r.length()==m.length,"Selected FASTA must preserve model order, safe output token and length");});
		require(n[0]==selected.size(),"Missing selected FASTA model");
	}
	private void loadMatrix(){
		final String path=required("matrix");final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{RrnaResourceIO.header(in,"#MITO_COVERAGE_MATRIX_V1\tSHORTER_QUERY_QUANTUM\tfull",path);
			for(String key:BINDINGS){RrnaResourceIO.header(in,"#input\t"+key+'\t'+bound(key),path);}
			RrnaResourceIO.header(in,MATRIX_HEADER,path);
			for(Unique u:uniques){
				for(Model m:models){final byte[] line=in.nextLine();require(line!=null,"Incomplete full coverage matrix");p.set(line);
					require(p.terms()==8 && p.termEquals(u.id,0) && p.termEquals(m.id,1) && p.parseInt(2)==u.length && p.parseInt(3)==m.length
						&& p.termEquals(u.length<=m.length?"MEMBER":"MODEL",4),"Matrix pair order, IDs, lengths and shorter-query direction must match exact input bindings");
					m.identity=p.parseFloat(5);m.trace=p.parseFloat(6);m.paired=p.parseInt(7);matrixRows++;
					require(Float.isFinite(m.identity) && m.identity>=0 && m.identity<=1 && Float.isFinite(m.trace) && m.trace>=0 && m.trace<=1 && m.paired>=0 && m.paired<=Math.min(u.length,m.length),"Invalid saved alignment evidence");
					require(m.original==m || m.identity==m.original.identity && m.trace==m.original.trace && m.paired==m.original.paired,"Exact model aliases must retain identical alignment evidence");
					if(m.step>0 && m.identity>=cut){m.coveredUnique++;m.coveredCopies+=u.copies;m.min=Math.min(m.min,u.length);m.max=Math.max(m.max,u.length);
						if(maxWindow>0){if(m.windowCovered!=null){m.windowCovered.set(u.index);u.ineligibleCovered=true;}else{u.feasibleCovered=true;}}}
				}
				Model first=null,best=null;for(Model m:selected){if(m.identity>=cut && first==null){first=m;}if(best==null || m.identity>best.identity){best=m;}}
				require(first==u.assigned && best==u.best && best.identity==u.bestIdentity,"Saved member assignment and best selected score must reproduce from the bound matrix");
			}
			require(in.nextLine()==null && matrixRows==(long)uniques.size()*models.size(),"Full matrix must contain every expected pair exactly once");
		}finally{require(!in.close(),"Matrix read failed");}
		for(Model m:selected){require(m.coveredUnique>=m.newUnique && m.coveredCopies>=m.newCopies && m.min>0 && m.max>=m.min,"Full coverage set must include every newly assigned member");m.lower=lowerBound(m.min);m.upper=upperBound(m.max);}
	}
	/** Exact rational floor/ceiling avoids binary floating-point rounding at integers. */
	public static int lowerBound(int min){require(min>0,"Observed minimum length must be positive");final long n=4L*min/5;require(n>0,"Derived minimum cannot be zero");return (int)n;}
	public static int upperBound(int max){require(max>0,"Observed maximum length must be positive");final long n=(6L*max+4)/5;require(n<=Integer.MAX_VALUE,"Derived maximum exceeds supported genomic coordinates");return (int)n;}
	private void write(){
		final ByteStreamWriter bounds=RrnaResourceIO.writer(prefix+".bounds.tsv"),profile=RrnaResourceIO.writer(prefix+".profile.tsv"),summary=RrnaResourceIO.writer(prefix+".summary.tsv"),inputs=RrnaResourceIO.writer(prefix+".inputs.tsv");final ByteBuilder b=new ByteBuilder();
		try{bounds.println(BOUNDS_HEADER);profile.println("model\tsource_model\tmodel_length\tcovered_copies\tcovered_unique\tmember_min\tmember_max\tmin\tmax\tcut\tscope");
			for(Model m:selected){bounds.print(b.clear().append(m.output).tab().append(m.length).tab().append(m.lower).tab().append(m.upper).nl());profile.print(b.clear().append(m.output).tab().append(m.id).tab().append(m.length).tab().append(m.coveredCopies).tab().append(m.coveredUnique).tab().append(m.min).tab().append(m.max).tab().append(m.lower).tab().append(m.upper).tab().appendSlow(cut).append("\tALL_COVERED_MEMBERS\n"));}
			summary.println("selected_models\tcopies\tunique_sequences\tmatrix_rows\tcut\tnew_alignment_calls\tscope\tinterpretation");summary.print(b.clear().append(selected.size()).tab().append(copies).tab().append(uniques.size()).tab().append(matrixRows).tab().appendSlow(cut).append("\t0\tALL_COVERED_MEMBERS\tSHORTER_QUERY_COVERAGE_NOT_CALLER_RECALL\n"));
			inputs.println("role\tpath\tsha80");for(String key:BINDINGS){input(inputs,b,key);}for(String key:INPUTS){input(inputs,b,key);}
		}finally{RrnaResourceIO.finish(bounds,profile,summary,inputs);}
		if(maxWindow>0){writeWindows();}
		System.err.println("MITO_LENGTH_PROFILE_PASS selected="+selected.size()+" copies="+copies+" unique="+uniques.size()+" matrix_rows="+matrixRows+" cut="+cut+" new_alignments=0");
	}
	/** Identity-floor feasibility only: actual seed/window/rank/length gates and
	 * original16S-model coverage are not predicted by these saved shorter-query fits. */
	private void writeWindows(){
		require(copyUnique.size==copies && copyLineage.size==copyUnique.size,"Window lineages must retain every original member copy");
		final ByteStreamWriter summary=RrnaResourceIO.writer(prefix+".window_summary.tsv"),lineages=RrnaResourceIO.writer(prefix+".window_lineages.tsv");final ByteBuilder b=new ByteBuilder();
		try{summary.println("max_window\tcut\tscope\tmodel\tmodel_length\tidentity_min_length\tselected_models\tineligible_models\tcovered_copies\tcovered_unique\tpopulation_copies\tpopulation_unique\tinterpretation");
			lineages.println("max_window\tcut\tscope\tmodel\tphylum\tcovered_copies\tcovered_unique\tunique_counts_additive_across_lineages\tinterpretation");
			for(Model m:selected){if(m.windowCovered!=null){windowScope(summary,lineages,b,"INELIGIBLE_MODEL",m,false);}}
			windowScope(summary,lineages,b,"INELIGIBLE_UNION",null,false);
			windowScope(summary,lineages,b,"INELIGIBLE_ONLY",null,true);
		}finally{RrnaResourceIO.finish(summary,lineages);}
		System.err.println("MITO_WINDOW_PROFILE_PASS maxwindow="+maxWindow+" selected="+selected.size()+" ineligible_models="+ineligibleModels+" new_alignments=0");
	}
	private void windowScope(ByteStreamWriter summary,ByteStreamWriter lineages,ByteBuilder b,String scope,Model model,boolean only){
		assert(maxWindow>0 && (model==null || model.windowCovered!=null)) : "Window scopes require an explicit limit and an ineligible model or union";
		final WindowGroup total=new WindowGroup();final WindowGroup[] groups=new WindowGroup[lineageNames.size()];
		for(int i=0;i<copyUnique.size;i++){final Unique u=uniques.get(copyUnique.get(i));
			final boolean keep=model==null ? u.ineligibleCovered && (!only || !u.feasibleCovered) : model.windowCovered.get(u.index);
			if(!keep){continue;}total.add(u.index);final int lineage=copyLineage.get(i);if(groups[lineage]==null){groups[lineage]=new WindowGroup();}groups[lineage].add(u.index);
		}
		require(model==null || total.copies==model.coveredCopies && total.unique.cardinality()==model.coveredUnique,"Per-model window coverage must equal its full covered set, never only greedy gains");
		long lineageCopies=0;final BitSet union=new BitSet();for(WindowGroup g:groups){if(g!=null){lineageCopies+=g.copies;union.or(g.unique);}}
		require(lineageCopies==total.copies && union.equals(total.unique),"Lineage copies and union of lineage unique IDs must conserve the global scope");
		final String name=model==null ? "ALL" : model.output;
		b.clear().append(maxWindow).tab().appendSlow(cut).tab().append(scope).tab().append(name).tab();
		if(model==null){b.append("NA\tNA");}else{b.append(model.length).tab().append(model.identityMinimum);}
		summary.print(b.tab().append(selected.size()).tab().append(ineligibleModels).tab().append(total.copies).tab().append(total.unique.cardinality()).tab().append(copies).tab().append(uniques.size()).tab().append(WINDOW_INTERPRETATION).nl());
		for(int i=0;i<groups.length;i++){final WindowGroup g=groups[i];if(g==null){continue;}
			lineages.print(b.clear().append(maxWindow).tab().appendSlow(cut).tab().append(scope).tab().append(name).tab().append(lineageNames.get(i)).tab().append(g.copies).tab().append(g.unique.cardinality()).append("\tfalse\t").append(WINDOW_INTERPRETATION).nl());}
	}
	private void input(ByteStreamWriter out,ByteBuilder b,String key){out.print(b.clear().append(key).tab().append(required(key)).tab().append(bound(key)).nl());}
	private interface Row {void accept(LineParser1 p);}
	private static void read(String path,String header,int terms,Row action){final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');try{RrnaResourceIO.header(in,header,path);for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==terms,"Unexpected table schema: "+path);action.accept(p);}}finally{require(!in.close(),"Table read failed: "+path);}}
	private static boolean bit(LineParser1 p,int column){require(p.termEquals("true",column) || p.termEquals("false",column),"Member partial flag must be explicit");return p.termEquals("true",column);}
	private String required(String key){final String v=flags.get(key);require(v!=null && !v.isEmpty(),"Missing argument: "+key);return v;}
	private String value(String key,String fallback){final String v=flags.get(key);return v==null?fallback:v;}
	private String bound(String key){final String path=required(key),hash=pins.get(path);require(hash!=null,"Missing producer-verified input: "+key+"="+path);return hash;}
	private static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final class Unique{Unique(String i,int n,int x){id=i;length=n;index=x;}final String id;final int length,index;long copies;Model assigned,best;float bestIdentity;boolean ineligibleCovered,feasibleCovered;}
	private static final class WindowGroup{long copies;final BitSet unique=new BitSet();void add(int index){assert(index>=0):"Unique indexes come from the validated ordered FASTA";copies++;unique.set(index);}}
	private static final class Model{Model(LineParser1 p){id=p.parseString(0);output=p.parseString(1);length=p.parseInt(2);canonical=p.parseString(3);require(length>0 && output.matches("mitochondria_coverage_[A-Za-z0-9_]+"),"Model must have a positive length and loader-safe output name");}
		final String id,output,canonical;final int length;Model original;int step,newUnique,assignedUnique,coveredUnique,min=Integer.MAX_VALUE,max,lower,upper,paired,identityMinimum;long newCopies,assignedCopies,coveredCopies;float identity,trace;BitSet windowCovered;}
	public static final String MATRIX_HEADER="unique_id\tmodel_id\tmember_length\tmodel_length\tquery_side\tidentity\ttrace_identity\tpaired_columns";
	public static final String BOUNDS_HEADER="model\tmodel_length\tmin_length\tmax_length";
	public static final String[] BINDINGS={"corpus","ssumodels","nohitmodels","ssumembers","nohitmembers","panelcases","panelgenes","truth","tree"};
	private static final String[] INPUTS={"matrix","unique","models","members","selection","summary","selected"};
	private static final String KNOWN="|corpus|ssumodels|nohitmodels|ssumembers|nohitmembers|panelcases|panelgenes|truth|tree|matrix|unique|models|members|selection|summary|selected|verified|outprefix|cut|expectedselected|maxwindow|";
	private static final String WINDOW_INTERPRETATION="IDENTITY_FLOOR_ONLY_RESOURCE_COVERAGE_NOT_CALLER_FN";
	private final ObjectMap<String,String> flags=new ObjectMap<String,String>(String.class,String.class);private ObjectMap<String,String> pins;
	private final ArrayList<Model> models=new ArrayList<Model>(),selected=new ArrayList<Model>();private final ObjectMap<String,Model> byModel=new ObjectMap<String,Model>(String.class,Model.class);
	private final ArrayList<Unique> uniques=new ArrayList<Unique>();private final ObjectMap<String,Unique> byUnique=new ObjectMap<String,Unique>(String.class,Unique.class);
	private final ObjectIntMap<String> lineageIndex=new ObjectIntMap<String>(String.class);private final ArrayList<String> lineageNames=new ArrayList<String>();
	private final IntList copyUnique=new IntList(),copyLineage=new IntList();
	private String prefix;private float cut;private int expectedSelected,canonicalModels,maxWindow=-1,ineligibleModels;private long selectionUniverse=-1,copies,partialCopies,panelCopies,matrixRows;
}
