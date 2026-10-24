package prok;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectIntMap;
import map.ObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;
import structures.IntHashMap;
import tax.TaxTree;

/** Single-linkage summaries of saved cluster identities, without alignments.
 * Every declared cut is reported; none is selected. Dominant source labels use
 * original-copy counts, preserve modal ties and do not establish RNA ancestry.
 * @author Ganyu */
public final class NcrnaMitoClusterMerge {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args,NcrnaMitoClusterMerge.class,false);
		try{final NcrnaMitoClusterMerge x=new NcrnaMitoClusterMerge();x.parse(pp.args);x.process();
			pp.outstream.println("MITO_CLUSTER_MERGE_PASS arm="+x.arm+" members="+x.expected+" models="+x.expectedModels+" cuts="+x.cuts.length+" unassigned="+x.unassigned+" linkage=single no_alignments=true no_cut_selected=true");
		}finally{Shared.closeStream(pp.outstream);}
	}
	private void parse(String[] args){
		final Parser p=new Parser();String verified=null;
		for(String arg:args){final String[] kv=Parse.splitOnFirst(arg,'=');final String k=kv[0].toLowerCase(Locale.ROOT),v=kv.length>1?kv[1]:null;
			if(k.equals("pairs")){pairs=v;}else if(k.equals("tree")){treeFile=v;}else if(k.equals("arm")){arm=v;}
			else if(k.equals("expected")){expected=Parse.parseIntKMG(v);}else if(k.equals("expectedmodels")){expectedModels=Parse.parseIntKMG(v);}
			else if(k.equals("verified")){verified=v;}else if(k.equals("cuts")){cuts=parseCuts(v);}
			else if(k.equals("outtree")){treeOutput=v;}
			else if(!p.parse(arg,k,v)){throw new IllegalArgumentException("Unknown flag: "+arg);}
		}
		membersFile=p.in1;output=p.out1;require(arm!=null && expected>0 && expectedModels>0 && expectedModels<=46340 && !p.append && p.maxReads<0,"Require explicit arm/denominators and a complete matrix; no append/read limit");
		require(verified!=null && output!=null && treeOutput!=null && !output.equals(treeOutput),"Require verified input manifest and distinct fresh out=table/outtree=forest");final ObjectMap<String,String> pins=EukPanelGffValidation.hashes(verified);
		for(String path:new String[]{membersFile,pairs,treeFile}){require(path!=null && pins.get(path)!=null,"Missing pinned saved input: "+path);}
		require(treeFile.indexOf(',')<0,"Use one explicitly pinned taxonomy tree, not an auxiliary path list");RrnaResourceIO.fresh(output,"");RrnaResourceIO.fresh(treeOutput,"");
	}
	private void process(){
		readPairs();readMembers();assignLineages();final String[] names=models.toArray(new String[models.size()]);final int[][] forest=maximumSpanningForest(identity,names);
		final ByteStreamWriter out=RrnaResourceIO.writer(output),tree=RrnaResourceIO.writer(treeOutput);final ByteBuilder b=new ByteBuilder();
		try{tree.println("arm\tmerge_index\tleft_model\tright_model\tidentity\tordering");
			for(int i=0;i<forest.length;i++){final int a=forest[i][0],c=forest[i][1];tree.print(b.clear().append(arm).tab().append(i+1).tab().append(names[a]).tab().append(names[c]).tab().append(identity[a][c],9).append("\tIDENTITY_DESC_THEN_LEXICAL_ENDPOINT_NAMES").nl());}
			out.println("arm\tidentity_cut\tmerged_groups\tlargest_rank\tgroup_key\tmodels\tmembers\tunique_sequences\tarm_members\tassigned_members\tunassigned_members\tunassigned_unique\tdominant_product\tproduct_support\tproduct_ties\tdominant_phylum\tphylum_support\tphylum_ties\tdominant_class\tclass_support\tclass_ties\tdominant_best_cm\tbest_cm_support\tbest_cm_ties\tlinkage\tlabel_weighting");
			for(double cut:cuts){final int[] roots=components(identity,cut);final Group[] byRoot=new Group[models.size()];final ArrayList<Group> groups=new ArrayList<Group>();
				for(int i=0;i<models.size();i++){final int root=roots[i];Group g=byRoot[root];if(g==null){g=new Group(models.get(i));byRoot[root]=g;groups.add(g);}g.models++;if(models.get(i).compareTo(g.key)<0){g.key=models.get(i);}}
				int assigned=0;for(Member m:members){if(m.model<0){continue;}byRoot[roots[m.model]].add(m);assigned++;}
				require(assigned+unassigned==expected,"Merged groups plus explicit unassigned members must conserve the arm");
				groups.sort(Comparator.comparingInt((Group g)->g.copies).reversed().thenComparing(Comparator.comparingInt((Group g)->g.unique.size()).reversed()).thenComparing(g->g.key));
				int total=0;for(Group g:groups){require(g.copies>0,"Every native model must retain assigned members");total+=g.copies;}require(total==assigned,"Component membership cannot duplicate or omit a source copy");
				for(int rank=0;rank<Math.min(3,groups.size());rank++){final Group g=groups.get(rank);
					b.clear().append(arm).tab().append(cut,2).tab().append(groups.size()).tab().append(rank+1).tab().append(g.key).tab().append(g.models).tab().append(g.copies).tab().append(g.unique.size()).tab().append(expected).tab().append(assigned).tab().append(unassigned).tab().append(unassignedUnique.size());
					g.product.append(b);g.phylum.append(b);g.clazz.append(b);g.cm.append(b);out.print(b.append("\tSINGLE\tORIGINAL_MEMBER_COPY").nl());
				}
			}
		}finally{RrnaResourceIO.finish(out,tree);}
	}
	private void readPairs(){
		identity=new double[expectedModels][expectedModels];final boolean[][] seen=new boolean[expectedModels][expectedModels];
		final ByteFile in=RrnaResourceIO.open(pairs);final LineParser1 p=new LineParser1('\t');long count=0;
		try{RrnaResourceIO.header(in,PAIR_HEADER,pairs);
			for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==14 && p.termEquals(arm,0),"Pair rows must have the declared schema and arm");
				if(p.termEquals("UNIVERSAL",1)){continue;}require(p.termEquals("CLUSTER_PAIR",1),"Unknown saved pair kind");
				final int a=model(p.parseString(2)),b=model(p.parseString(3));final double value=p.parseDouble(7);
				require(Double.isFinite(value) && value>=0 && value<=1 && !seen[a][b],"Pair identities must be valid and unordered pairs unique");
				seen[a][b]=seen[b][a]=true;identity[a][b]=identity[b][a]=value;count++;
			}
		}finally{require(!in.close(),"Saved matrix read failed");}
		require(models.size()==expectedModels && count==(long)expectedModels*(expectedModels+1)/2,"Saved matrix must contain every unordered pair including every diagonal");
		for(int i=0;i<expectedModels;i++){for(int j=i;j<expectedModels;j++){require(seen[i][j],"Missing pair cannot be treated as low identity");}}
	}
	private int model(String name){
		require(name!=null && !name.isEmpty(),"Matrix model names must be present");int index=modelIds.get(name);
		if(index<0){require(models.size()<expectedModels,"More model names than declared");final String token=RrnaResourceIO.first(name);require(modelTokens.get(token)<0,"Duplicate first-token model ID");index=models.size();models.add(name);modelIds.put(name,index);modelTokens.put(token,index);}return index;
	}
	private void readMembers(){
		final ByteFile in=RrnaResourceIO.open(membersFile);final LineParser1 p=new LineParser1('\t');final ObjectIntMap<String> ids=new ObjectIntMap<String>(String.class),unique=new ObjectIntMap<String>(String.class);final int[] counts=new int[expectedModels];
		try{RrnaResourceIO.header(in,MEMBER_HEADER,membersFile);
			for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==22 && p.termEquals(arm,0),"Saved member rows must match the declared schema/arm");final String id=p.parseString(1);require(!id.isEmpty() && ids.get(id)<0,"Missing or duplicate source locus: "+id);ids.put(id,1);
				final Member m=new Member();final String uid=p.parseString(2);require(!uid.isEmpty(),"Every source copy needs its retained unique-sequence ID");m.unique=unique.get(uid);if(m.unique<0){m.unique=unique.size();unique.put(uid,m.unique);}
				m.taxid=p.parseInt(3);require(m.taxid>0,"Source host taxid must remain positive");m.product=p.parseString(7);m.cm=p.parseString(12);
				if(p.termEquals("UNASSIGNED",16)){require((p.termEquals("small",14) || p.termEquals("orphan",14)) && p.termEquals("NA",17),"Unassigned rows require a recognized native status and no invented identity");unassigned++;unassignedUnique.put(m.unique,1);}
				else{require(p.termEquals("kept",14),"Assigned members must retain native kept status");m.model=modelIds.get(p.parseString(16));require(m.model>=0,"Assigned model absent from pair matrix");counts[m.model]++;}
				members.add(m);
			}
		}finally{require(!in.close(),"Saved member read failed");}
		require(members.size()==expected,"Merge summaries must retain the complete arm denominator");
		for(int i=0;i<models.size();i++){require(counts[i]>0 && counts[i]==Integer.parseInt(RrnaResourceIO.field(models.get(i),"n")),"Native model n= differs from assigned member count: "+models.get(i));}
	}
	private void assignLineages(){
		assert(!members.isEmpty()):"Source membership must be validated before taxonomy grouping";final TaxTree tree;final boolean old=TaxTree.LOAD_FAKE_TREE;
		try{TaxTree.LOAD_FAKE_TREE=false;tree=TaxTree.loadTaxTreeFile(treeFile);}finally{TaxTree.LOAD_FAKE_TREE=old;}
		final IntHashMap cache=new IntHashMap();final ArrayList<String[]> lineages=new ArrayList<String[]>();
		for(Member m:members){int index=cache.get(m.taxid);if(index<0){index=lineages.size();lineages.add(new String[]{NcrnaMitoClusterEvaluation.lineage(tree,m.taxid,TaxTree.PHYLUM_E),NcrnaMitoClusterEvaluation.lineage(tree,m.taxid,TaxTree.CLASS_E)});cache.put(m.taxid,index);}m.phylum=lineages.get(index)[0];m.clazz=lineages.get(index)[1];}
	}
	/** Each edge at or above the cut joins components; bridges are deliberately transitive. */
	public static int[] components(double[][] matrix,double cut){
		require(matrix!=null && matrix.length>0 && Double.isFinite(cut) && cut>=0 && cut<=1,"Single linkage requires a nonempty square matrix and explicit valid cut");final int[] parent=new int[matrix.length];for(int i=0;i<parent.length;i++){require(matrix[i]!=null && matrix[i].length==parent.length,"Identity matrix must be square");parent[i]=i;}
		for(int i=0;i<parent.length;i++){for(int j=i+1;j<parent.length;j++){require(Double.isFinite(matrix[i][j]) && matrix[i][j]>=0 && matrix[i][j]<=1 && matrix[i][j]==matrix[j][i],"Single linkage consumes a complete symmetric saved matrix");if(matrix[i][j]>=cut){final int a=find(parent,i),b=find(parent,j);if(a!=b){parent[Math.max(a,b)]=Math.min(a,b);}}}}
		for(int i=0;i<parent.length;i++){parent[i]=find(parent,i);}return parent;
	}
	/** Kruskal maximum-spanning forest; only actual union edges are emitted.
	 * Equal identities use lexical, canonically ordered endpoint names, independent
	 * of input model order. Retaining edges at or above any cut reproduces single linkage. */
	public static int[][] maximumSpanningForest(double[][] matrix,String[] names){
		require(matrix!=null && names!=null && matrix.length==names.length && names.length>0,"Forest requires the exact nonempty named matrix");
		final int n=names.length;final int[] parent=new int[n];final ObjectIntMap<String> used=new ObjectIntMap<String>(String.class);final ArrayList<Edge> edges=new ArrayList<Edge>();
		for(int i=0;i<n;i++){require(matrix[i]!=null && matrix[i].length==n && names[i]!=null && !names[i].isEmpty() && used.get(names[i])<0,"Forest model names must be unique and matrix square");used.put(names[i],i);parent[i]=i;}
		for(int i=0;i<n;i++){for(int j=i+1;j<n;j++){require(Double.isFinite(matrix[i][j]) && matrix[i][j]>=0 && matrix[i][j]<=1 && matrix[i][j]==matrix[j][i],"Forest requires valid symmetric saved identities");edges.add(names[i].compareTo(names[j])<0?new Edge(i,j,matrix[i][j]):new Edge(j,i,matrix[i][j]));}}
		edges.sort(Comparator.comparingDouble((Edge e)->e.identity).reversed().thenComparing(e->names[e.a]).thenComparing(e->names[e.b]));
		final int[][] result=new int[n-1][2];int count=0;
		for(Edge e:edges){final int a=find(parent,e.a),b=find(parent,e.b);if(a==b){continue;}parent[Math.max(a,b)]=Math.min(a,b);result[count][0]=e.a;result[count][1]=e.b;count++;if(count==n-1){break;}}
		require(count==n-1,"A complete finite pair matrix must span every model");return result;
	}
	private static int find(int[] parent,int i){assert(i>=0 && i<parent.length):"Component index comes from the checked model roster";while(parent[i]!=i){parent[i]=parent[parent[i]];i=parent[i];}return i;}
	private static double[] parseCuts(String text){require(text!=null && !text.isEmpty(),"Cut list must be explicit");final String[] fields=text.split(",");final double[] values=new double[fields.length];for(int i=0;i<values.length;i++){values[i]=Double.parseDouble(fields[i]);require(Double.isFinite(values[i]) && values[i]>=0 && values[i]<=1 && (i==0 || values[i]>values[i-1]),"Cuts must be finite, ascending, distinct fractions");}return values;}
	private static void require(boolean value,String message){if(!value){throw new IllegalArgumentException(message);}}
	private static final class Member{int model=-1,unique,taxid;String product,phylum,clazz,cm;}
	private static final class Edge{Edge(int x,int y,double d){a=x;b=y;identity=d;}final int a,b;final double identity;}
	private static final class Group{
		Group(String k){key=k;}void add(Member m){assert(m.model>=0):"Only assigned copies belong to a consensus component";copies++;unique.put(m.unique,1);product.add(m.product);phylum.add(m.phylum);clazz.add(m.clazz);cm.add(m.cm);}
		String key;int models,copies;final IntHashMap unique=new IntHashMap();final Mode product=new Mode(),phylum=new Mode(),clazz=new Mode(),cm=new Mode();
	}
	private static final class Mode{
		void add(String label){assert(label!=null):"Dominance preserves an explicit source label or taxonomy-gap marker";final int n=counts.get(label);if(n<0){labels.add(label);counts.put(label,1);}else{counts.put(label,n+1);}}
		void append(ByteBuilder b){require(!labels.isEmpty(),"A retained group needs source-label observations");int max=0;final ArrayList<String> tied=new ArrayList<String>();for(String label:labels){final int n=counts.get(label);if(n>max){max=n;tied.clear();}if(n==max){tied.add(label);}}tied.sort(String::compareTo);b.tab().append(String.join(" | ",tied)).tab().append(max).tab().append(tied.size());}
		final ObjectIntMap<String> counts=new ObjectIntMap<String>(String.class);final ArrayList<String> labels=new ArrayList<String>();
	}
	private String membersFile,pairs,treeFile,output,treeOutput,arm;private int expected,expectedModels,unassigned;private double[][] identity;
	private double[] cuts={.60,.65,.70,.75,.80,.85,.90};private final ArrayList<String> models=new ArrayList<String>();private final ArrayList<Member> members=new ArrayList<Member>();
	private final ObjectIntMap<String> modelIds=new ObjectIntMap<String>(String.class),modelTokens=new ObjectIntMap<String>(String.class);private final IntHashMap unassignedUnique=new IntHashMap();
	public static final String PAIR_HEADER="arm\tkind\tleft_id\tright_id\tleft_length\tright_length\tquery_side\tquantum_identity\ttrace_identity\tref_start0\tref_stop0\tpaired_columns\tleft_paired_fraction\tright_paired_fraction";
	public static final String MEMBER_HEADER="arm\tlocus_id\tunique_id\tsource_taxid\tsource_length\tpartial\tpseudo\tsource_product\tlength_bin\tlineage_rank\tlineage\tuniversal_preference\tbest_ssu_cm\tcm_concordance\tassignment_status\tcluster_key\tassigned_model\tquantum_identity\ttrace_identity\tpaired_columns\tmember_paired_fraction\tmodel_paired_fraction";
}
