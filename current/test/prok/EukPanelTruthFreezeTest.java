package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import prok.EukPanelTruthFreeze;
import structures.ByteBuilder;

/** Freeze must preserve all bytes outside the status field and reject wrong inputs.
 * @author Raiden */
public final class EukPanelTruthFreezeTest {
	public static void main(String[] args){
		require(args.length==1 && args[0].startsWith("outdir="),"Fresh outdir required");final String root=args[0].substring(7);require(new File(root).mkdir(),"Fresh fixture root required");
		final String[] p=new String[35];java.util.Arrays.fill(p,"x");p[0]="17";p[32]="FREEZE_CANDIDATE_NOT_FROZEN";p[34]="1.5_WITH_HOMOLOG_RULE_WITHDRAWN";
		final String input=EukPanelTruthFreeze.HEADER+'\n'+String.join("\t",p)+'\n';write(root+"/input.tsv",input);run(root,"input.tsv","truth.tsv.gz",1);
		require(read(root+"/truth.tsv.gz").equals(input.replace("FREEZE_CANDIDATE_NOT_FROZEN","TRUTH_V1")),"Only status field may change");
		write(root+"/duplicate.tsv",input+String.join("\t",p)+'\n');fail(root,"duplicate.tsv","duplicate.out",2,"Duplicate locus");
		fail(root,"input.tsv","wrongcount.out",2,"denominator");fail(root,"truth.tsv.gz","refreeze.out",1,"Only the accepted");
		final String v2=input.replace("FREEZE_CANDIDATE_NOT_FROZEN",EukPanelTruthFreeze.V2_CANDIDATE).replace("1.5_WITH_HOMOLOG_RULE_WITHDRAWN",EukPanelTruthFreeze.V2_POLICY);
		write(root+"/v2candidate.tsv",v2);
		EukPanelTruthFreeze.main(new String[]{"in="+root+"/v2candidate.tsv","out="+root+"/v2truth.tsv.gz","expectedloci=1","version=2","policy="+EukPanelTruthFreeze.V2_POLICY});
		require(read(root+"/v2truth.tsv.gz").equals(v2.replace(EukPanelTruthFreeze.V2_CANDIDATE,"TRUTH_V2")),"Version2 freeze changes only the status field");
		fail(root,"v2candidate.tsv","wrongversion.out",1,"Only the accepted");
		boolean badPolicy=false;try{EukPanelTruthFreeze.main(new String[]{"in="+root+"/input.tsv","out="+root+"/badpolicy.out","expectedloci=1","version=2"});}
		catch(IllegalArgumentException e){badPolicy=e.getMessage().contains("policy/version");}require(badPolicy,"Version2 must explicitly select its reviewed policy");
		System.out.println("PASS status-only byte preservation, unique IDs, exact denominator and rejection of refreezing");
	}
	static void run(String root,String in,String out,int n){EukPanelTruthFreeze.main(new String[]{"in="+root+'/'+in,"out="+root+'/'+out,"expectedloci="+n});}
	static void fail(String root,String in,String out,int n,String reason){boolean failed=false;try{run(root,in,out,n);}catch(IllegalArgumentException e){require(e.getMessage().contains(reason),"Wrong failure: "+e);failed=true;}require(failed,"Expected rejection");}
	static void write(String path,String text){final ByteStreamWriter w=new ByteStreamWriter(path,false,false,true);w.start();w.print(text);require(!w.poisonAndWait(),"Fixture write failed");}
	static String read(String path){final ByteFile r=ByteFile.makeByteFile(FileFormat.testInput(path,FileFormat.TEXT,null,true,true));final ByteBuilder b=new ByteBuilder();try{for(byte[] line;(line=r.nextLine())!=null;){b.append(line).nl();}}finally{require(!r.close(),"Fixture read failed");}return new String(b.toBytes(),StandardCharsets.UTF_8);}
	static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
}
