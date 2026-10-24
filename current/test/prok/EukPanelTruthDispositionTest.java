package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.LineParser1;
import prok.EukPanelTruthDisposition;
import structures.ByteBuilder;

/** Cross-package fixtures exercise actual exported membership and metric eligibility.
 * @author Raiden */
public final class EukPanelTruthDispositionTest {
	public static void main(String[] args){
		require(args.length==1 && args[0].startsWith("outdir="),"Expected fresh outdir");final String root=args[0].substring(7);require(new File(root).mkdir(),"Fresh fixture directory required");
		final int[] starts={100,100,300,150,600,590,690},ends={399,199,399,249,699,639,790};
		final String[] models={"RF01960","RF00177","RF01959","RF02543","RF02541","RF02540","RF02543"};
		final String[] families={"SSU","SSU","SSU","LSU","LSU","LSU","LSU"};
		final int[] bits={300,100,90,100,110,70,80};
		final ByteBuilder loci=new ByteBuilder(),members=new ByteBuilder();loci.append(EukPanelTruthDisposition.INPUT_HEADER).nl();members.append(EukPanelTruthDisposition.MEMBER_HEADER).nl();
		for(int i=0;i<7;i++){final int id=i+1;final String flags=id==1 || id==5?"SPANNING_AMBIGUOUS":id==2 || id==4?"AMBIGUOUS_CROSSFOLD":"NONE";
			loci.append(locus(id,starts[i],ends[i],families[i],models[i],bits[i],flags,id==2));members.append(member(id,starts[i],ends[i],families[i],models[i],bits[i]));}
		write(root+"/candidate.tsv",loci.toString());write(root+"/members.tsv",members.toString());
		final String cases=EukPanelTruthDisposition.CASE_HEADER+"\n1\tSPANNING_SEED\tA\tchr\t+\tSSU\trun\t1\t100\t399\trun\t2\t100\t199\trun\t3\t300\t399\n"
			+"5\tSPANNING_SEED\tA\tchr\t+\tLSU\trun\t5\t600\t699\trun\t6\t590\t639\trun\t7\t690\t790\n";
		write(root+"/cases.tsv",cases);run(root,"members.tsv","pass");
		require(read(root+"/pass.summary.tsv").contains("7\t5\t2\t7\t1\t2\t2\tFREEZE_CANDIDATE_NOT_FROZEN"),"Exactly two whole witness loci absorbed; seven members conserved; two both-end exclusions");
		final ByteFile in=open(root+"/pass.truth_candidate.tsv.gz");final LineParser1 p=new LineParser1('\t');int n=0;
		try{in.nextLine();for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==35,"Disposition output schema");final int id=p.parseInt(0);
			if(id==1){require(p.parseLong(4)==100 && p.parseLong(5)==399 && p.parseInt(18)==3 && p.termEquals("RF00177;RF01959;RF01960",19) && p.termEquals("RF01960",25),"Keep euk seed span; recompute member set and subtype");
				require(p.parseString(28).contains("SPANNING_ENDS_PATTERN") && p.parseString(28).contains("AMBIGUOUS_CROSSFOLD"),"Ends disposition and absorbed crossfold uncertainty survive");
				require(p.termEquals("0",29) && p.termEquals("0",30),"Do not OR historical evidence beyond the surviving seed span");}
			if(id==5){require(p.parseLong(4)==600 && p.parseLong(5)==699 && p.parseInt(18)==1 && p.termEquals("SPANNING_AMBIGUOUS",28),"Exception retains its own locus and members");}
			require(p.termEquals("1",31) && p.termEquals(id==1 || id==4?"0":"1",33),"Crossfold excludes only BOTH ends, not exact5/exact3 eligibility");n++;
		}}finally{require(!in.close(),"Read candidate fixture");}require(n==5,"No unexpected locus loss");
		final ByteFile mi=open(root+"/pass.locus_members.tsv.gz");int row=0;try{mi.nextLine();for(byte[] line;(line=mi.nextLine())!=null;){p.set(line);row++;require(p.terms()==43 && p.parseLong(4)==row,"Original hit order and identities conserved");
			require(p.parseInt(40)==(row<=3?1:row),"Every hit follows the explicit old-to-new map");require(p.termEquals(row==2 || row==3?"ENDS_WITNESS":"SEED",41),"One surviving SEED and retained witness roles");}}
		finally{require(!mi.close(),"Read members fixture");}require(row==7,"Seven raw members preserved");
		require(read(root+"/pass.merge_map.tsv").split("\n").length==3,"Two explicit absorbed-locus mappings");
		write(root+"/duplicate.tsv",members.toString()+member(7,starts[6],ends[6],families[6],models[6],bits[6]));
		boolean failed=false;try{run(root,"duplicate.tsv","fail");}catch(IllegalArgumentException e){require(e.getMessage().contains("Repeated original member"),"Expected identity failure, got "+e);failed=true;}require(failed,"Duplicate raw members must fail before output");
		System.out.println("PASS witness-locus absorption, conserved hit identities, fixed seed geometry, unchanged exception, historical-support scope and distinct boundary denominators");
	}
	static String locus(int id,int start,int end,String family,String accession,int bits,String flags,boolean hist){
		return id+"\tA\tchr\t+\t"+start+'\t'+end+'\t'+family+"\ttest\tNUCLEAR\t"+accession+"\trun\t/raw.tbl\t"+id+'\t'+bits+"\t1\t1\t0\tINTERNAL\t1\t"+accession+"\t0\t0\t"+flags+"\tDESCRIPTIVE_ONLY\t1\t"+accession+"\tADMITTED\tFULLSPAN_INTACT\t"+flags+'\t'+(hist?1:0)+'\t'+(hist?1:0)+"\t1\tFREEZE_CANDIDATE_NOT_FROZEN\n";
	}
	static String member(int id,int start,int end,String family,String accession,int bits){
		return "A\trun\t/raw.tbl\t00000000000000000000\t"+id+"\t00000000000000000000\t"+accession+'\t'+family+"\ttest\tchr\tchr\tNUCLEAR\tEXACT\tVALID_INTERVAL\t1\t100\t100\t100\t1\t"+start+'\t'+end+"\t+\tno\t!\t"+bits+"\t1e-10\t50\t"+(bits-50)+"\t0\ttest\t1000\tINTERNAL\ttrue\t0\t0\ttrue\t0\t0\t0\tREVIEW_ONLY\t"+id+"\tSEED\t0\n";
	}
	static void run(String root,String members,String name){EukPanelTruthDisposition.main(new String[]{"candidate="+root+"/candidate.tsv","members="+root+'/'+members,"cases="+root+"/cases.tsv","exception=5","expectedends=1","expectedloci=7","expectedmembers=7","outprefix="+root+'/'+name});}
	static ByteFile open(String path){return ByteFile.makeByteFile(FileFormat.testInput(path,FileFormat.TEXT,null,true,true));}
	static String read(String path){final ByteFile in=open(path);final ByteBuilder b=new ByteBuilder();try{for(byte[] line;(line=in.nextLine())!=null;){b.append(line).nl();}}finally{require(!in.close(),"Fixture read failure");}return new String(b.toBytes(),StandardCharsets.UTF_8);}
	static void write(String path,String text){final ByteStreamWriter out=new ByteStreamWriter(path,false,false,true);out.start();out.print(text);require(!out.poisonAndWait(),"Fixture write failure");}
	static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
}
