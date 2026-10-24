package prok;

import java.util.ArrayList;
import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import map.ObjectIntMap;
import map.ObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import stream.Writer;
import stream.WriterFactory;
import structures.ByteBuilder;
import tax.TaxTree;

/** Exports one exact single-linkage component, never a dominant-taxon filter.
 * Full source/member/model rosters are checked before any output is written.
 * The retained native cluster seeds preserve original n= copy weighting.
 * @author Brian, Raiden */
public final class NcrnaMitoGroupExport {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args,NcrnaMitoGroupExport.class,false);
		try{final NcrnaMitoGroupExport x=new NcrnaMitoGroupExport();x.parse(pp.args);x.process();
			pp.outstream.println("MITO_GROUP_EXPORT_PASS group="+x.group+" input="+x.expected+" selected="+x.selected.size()+" models="+x.seeds.size()+" cut="+x.cut+" exact_component=true taxon_filter=false");
		}finally{Shared.closeStream(pp.outstream);}
	}
	private void parse(String[] args){
		final Parser p=new Parser();String verified=null;
		for(String arg:args){final String[] kv=Parse.splitOnFirst(arg,'=');final String k=kv[0].toLowerCase(Locale.ROOT),v=kv.length>1?kv[1]:null;
			if(k.equals("members")){memberFile=v;}else if(k.equals("pairs")){pairFile=v;}else if(k.equals("models")){modelFile=v;}else if(k.equals("arm")){arm=v;}
			else if(k.equals("target")){target=v;}else if(k.equals("group")){group=v;}else if(k.equals("outprefix")){prefix=v;}else if(k.equals("cut")){cut=Double.parseDouble(v);}
			else if(k.equals("expected")){expected=Parse.parseIntKMG(v);}else if(k.equals("expectedmodels")){expectedModels=Parse.parseIntKMG(v);}else if(k.equals("selected")){expectedSelected=Parse.parseIntKMG(v);}else if(k.equals("verified")){verified=v;}
			else if(!p.parse(arg,k,v)){throw new IllegalArgumentException("Unknown flag: "+arg);}
		}
		input=p.in1;require(expected>0 && expectedModels>0 && expectedModels<=46340 && expectedSelected>0 && expectedSelected<=expected && Double.isFinite(cut) && cut>=0 && cut<=1 && arm!=null && target!=null && group!=null && prefix!=null && !p.append && p.maxReads<0,"Require complete source/model/selected denominators, explicit component key/cut, fresh output, and no read limit");
		require(verified!=null,"All source tables and FASTAs require explicit producer pins");final ObjectMap<String,String> pins=EukPanelGffValidation.hashes(verified);for(String path:new String[]{input,memberFile,pairFile,modelFile}){require(path!=null && pins.get(path)!=null,"Missing verified input: "+path);}
		Read.TO_UPPER_CASE=true;Shared.TRIM_READ_DESCRIPTION=false;RrnaResourceIO.fresh(prefix,SUFFIXES);
	}
	private void process(){
		RrnaResourceIO.stream(modelFile,r->{final String id=RrnaResourceIO.first(r.id);require(r.bases!=null && r.length()>0 && tokens.get(id)<0,"Every native model needs unique token and nonempty bases");tokens.put(id,models.size());modelIds.put(r.id,models.size());models.add(new Read(r.bases.clone(),null,r.id,models.size()));});
		require(models.size()==expectedModels,"Model count differs from the complete saved matrix roster");
		final double[][] matrix=readPairs();final int[] components=NcrnaMitoClusterMerge.components(matrix,cut);final int targetIndex=tokens.get(target);require(targetIndex>=0,"Requested group key is absent from the model library: "+target);
		final int component=components[targetIndex];String key=null;final boolean[] included=new boolean[models.size()];
		for(int i=0;i<models.size();i++){if(components[i]==component){included[i]=true;seeds.add(models.get(i));if(key==null || models.get(i).id.compareTo(key)<0){key=models.get(i).id;}}}
		require(RrnaResourceIO.first(key).equals(target),"Use the canonical lexical component key from the accepted cut table, not an arbitrary member model");
		final ObjectMap<String,Member> roster=readMembers(included);int[] seen={0};
		RrnaResourceIO.stream(input,r->{final String id=RrnaResourceIO.first(r.id);final Member m=roster.get(id);require(m!=null && !m.seen && r.bases!=null && r.length()==m.length,"Source FASTA must conserve each accepted member exactly once: "+id);
			require(TaxTree.parseTaxID(r.id)==m.taxid && RrnaResourceIO.field(r.id,"unique").equals(m.uid) && flag(RrnaResourceIO.field(r.id,"partial"))==m.partial && flag(RrnaResourceIO.field(r.id,"pseudo"))==m.pseudo,"Source taxid/unique/flags differ from the accepted member table: "+id);
			m.seen=true;seen[0]++;if(m.included){selected.add(new Read(r.bases.clone(),null,r.id,r.numericID));selectedMembers.add(m);}
		});require(seen[0]==expected && selected.size()==expectedSelected,"Exact component export must match both full input and locked selected counts");
		final NcrnaMitoConsensus.Data data=new NcrnaMitoConsensus.Data(selected);int min=Integer.MAX_VALUE,partial=0,pseudo=0;long bases=0;
		for(Read r:selected){min=Math.min(min,r.length());bases+=r.length();}for(Member m:selectedMembers){partial+=m.partial?1:0;pseudo+=m.pseudo?1:0;}
		writeFasta(prefix+SUFFIXES[0],selected);writeFasta(prefix+SUFFIXES[1],seeds);
		final ByteStreamWriter members=RrnaResourceIO.writer(prefix+SUFFIXES[2]),summary=RrnaResourceIO.writer(prefix+SUFFIXES[3]);final ByteBuilder b=new ByteBuilder();
		try{members.println("group\tlocus_id\tunique_id\ttaxid\tpartial\tpseudo\tnative_model\tcomponent_key\tcut");for(Member m:selectedMembers){members.print(b.clear().append(group).tab().append(m.id).tab().append(m.uid).tab().append(m.taxid).tab().append(m.partial).tab().append(m.pseudo).tab().append(models.get(m.model).id).tab().append(key).tab().append(cut,2).nl());}
			summary.println("group\tarm\tcomponent_key\tcut\tinput_copies\tselected_copies\tselected_unique\tselected_models\tpartial\tpseudo\tmin_length\tmax_length\tbases\tlargest_seed\ttaxon_filter");
			summary.print(b.clear().append(group).tab().append(arm).tab().append(key).tab().append(cut,2).tab().append(expected).tab().append(selected.size()).tab().append(data.unique.size()).tab().append(seeds.size()).tab().append(partial).tab().append(pseudo).tab().append(min).tab().append(data.maxLength).tab().append(bases).tab().append(NcrnaMitoConsensus.largestSeed(seeds).id).append("\tfalse").nl());
		}finally{RrnaResourceIO.finish(members,summary);}
	}
	private double[][] readPairs(){
		final double[][] matrix=new double[models.size()][models.size()];final boolean[][] seen=new boolean[models.size()][models.size()];final ByteFile in=RrnaResourceIO.open(pairFile);final LineParser1 p=new LineParser1('\t');long count=0;
		try{RrnaResourceIO.header(in,NcrnaMitoClusterMerge.PAIR_HEADER,pairFile);for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==14 && p.termEquals(arm,0),"Pair schema/arm changed");if(p.termEquals("UNIVERSAL",1)){continue;}require(p.termEquals("CLUSTER_PAIR",1),"Unexpected pair kind");
			final int a=modelIds.get(p.parseString(2)),b=modelIds.get(p.parseString(3));final double identity=p.parseDouble(7);require(a>=0 && b>=0 && !seen[a][b] && Double.isFinite(identity) && identity>=0 && identity<=1,"Duplicate/unknown pair or nonfinite identity");seen[a][b]=seen[b][a]=true;matrix[a][b]=matrix[b][a]=identity;count++;
		}}finally{require(!in.close(),"Pair read failed");}
		require(count==(long)models.size()*(models.size()+1)/2,"Complete unordered matrix including diagonals is required");return matrix;
	}
	private ObjectMap<String,Member> readMembers(boolean[] included){
		final ObjectMap<String,Member> rows=new ObjectMap<String,Member>(String.class,Member.class);final ByteFile in=RrnaResourceIO.open(memberFile);final LineParser1 p=new LineParser1('\t');final int[] counts=new int[models.size()];
		try{RrnaResourceIO.header(in,NcrnaMitoClusterMerge.MEMBER_HEADER,memberFile);for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==22 && p.termEquals(arm,0),"Member schema/arm changed");
			final Member m=new Member();m.id=p.parseString(1);m.uid=p.parseString(2);m.taxid=p.parseInt(3);m.length=p.parseInt(4);m.partial=flag(p.parseString(5));m.pseudo=flag(p.parseString(6));require(!m.id.isEmpty() && !m.uid.isEmpty() && m.taxid>0 && m.length>0 && rows.get(m.id)==null,"Missing/duplicate member identity or invalid source metadata");
			if(p.termEquals("kept",14)){m.model=modelIds.get(p.parseString(16));require(m.model>=0,"Assigned model absent from source library");counts[m.model]++;m.included=included[m.model];}
			else{require((p.termEquals("small",14) || p.termEquals("orphan",14)) && p.termEquals("UNASSIGNED",16),"Unknown unassigned status");}rows.put(m.id,m);
		}}finally{require(!in.close(),"Member read failed");}
		require(rows.size()==expected,"Every source copy must have an accepted evaluation row");for(int i=0;i<models.size();i++){require(counts[i]>0 && counts[i]==Integer.parseInt(RrnaResourceIO.field(models.get(i).id,"n")),"Native n= must conserve original member counts for every source model");}return rows;
	}
	private static void writeFasta(String path,ArrayList<Read> records){assert(!records.isEmpty()):"An accepted component needs original members and native seeds";final Writer out=WriterFactory.makeWriter(FileFormat.testOutput(path,FileFormat.FASTA,null,true,false,false,true),null,1);out.start();out.add(records,0);require(!out.poisonAndWait() && out.readsWritten()==records.size(),"FASTA export must preserve the full selected roster");}
	private static boolean flag(String s){require(s.equals("true") || s.equals("false"),"Source flags must remain explicit booleans");return s.equals("true");}
	private static void require(boolean value,String message){if(!value){throw new IllegalArgumentException(message);}}
	private static final class Member{String id,uid;int taxid,length,model=-1;boolean partial,pseudo,included,seen;}
	private String input,memberFile,pairFile,modelFile,arm,target,group,prefix;private int expected,expectedModels,expectedSelected;private double cut=Double.NaN;
	private final ObjectIntMap<String> modelIds=new ObjectIntMap<String>(String.class),tokens=new ObjectIntMap<String>(String.class);private final ArrayList<Read> models=new ArrayList<Read>(),seeds=new ArrayList<Read>(),selected=new ArrayList<Read>();private final ArrayList<Member> selectedMembers=new ArrayList<Member>();
	private static final String[] SUFFIXES={".fa.gz",".seeds.fa.gz",".members.tsv",".summary.tsv"};
}
