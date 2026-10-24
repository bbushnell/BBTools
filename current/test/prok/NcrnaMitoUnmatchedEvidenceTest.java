package test.prok;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import map.ObjectMap;
import prok.EukPanelTruthFreeze;
import prok.NcrnaMitoUnmatchedEvidence;
import prok.NcrnaP1MissInventory;

/** Original-row binding, zero-overlap organelle arms and honest origin/identity fixtures.
 * @author Brian, Ganyu */
public final class NcrnaMitoUnmatchedEvidenceTest {
	public static void main(String[] args)throws Exception{
		check(args.length==1,"Supply a fresh fixture directory");final Path dir=Paths.get(args[0]);Files.createDirectory(dir);
		final String hash="00000000000000000001";
		final Path nuclear=write(dir,"nuclear.fa",">chr1 original nuclear scaffold\n"+repeat("ACGT",20)+"\n");
		final Path mito=write(dir,"mito.fa",">mito1 original mitochondrion header\n"+repeat("ACGT",10)+"\n");
		final String old=gff("chr1",1,5,"+","16S","mitochondria_insect"),added=gff("chr1",10,25,"-","16S","mitochondria_insect");
		final Path p2base=write(dir,"p2base.gff","##gff-version 3\n"+old+"\n");
		final Path p2a=write(dir,"p2a.gff","##gff-version 3\n"+old+"\n");
		final Path p2b=write(dir,"p2b.gff","##gff-version 3\n"+old+"\n"+added+"\n");
		final Path p3base=write(dir,"p3base.gff","##gff-version 3\n"+gff("mito1",15,25,"+","euk18S","euk18S_mitochondria")+"\n");
		final Path p3a=write(dir,"p3a.gff","##gff-version 3\nmito1\tfixture\tCDS\t25\t35\t10\t+\t0\tproduct:CDS\nmito1\tfixture\tregion\t1\t40\t.\t+\t.\tID=region1\n");
		final Path p3b=write(dir,"p3b.gff","##gff-version 3\nmito1\tfixture\ttRNA\t21\t25\t.\t-\t.\tproduct:tRNA\n");
		final ArrayList<String> audit=new ArrayList<String>();audit.add(NcrnaP1MissInventory.AUDIT_HEADER);
		audit.add(audit("baseline","P2","N","chr1","NUCLEAR","CALL",p2base+":2","16S",1,5,"+","APPARENT_FP"));
		audit.add(audit("B","P2","N","chr1","NUCLEAR","CALL",p2b+":2","16S",1,5,"+","APPARENT_FP"));
		audit.add(audit("B","P2","N","chr1","NUCLEAR","CALL",p2b+":3","16S",10,25,"-","APPARENT_FP"));
		for(String arm:new String[]{"baseline","A","B"}){audit.add(audit(arm,"P3","M","mito1","MITOCHONDRION","TRUTH","LOCUS:1","RF00177",20,30,"+","FN"));}
		final Path auditPath=dir.resolve("audit.tsv");Files.write(auditPath,audit,StandardCharsets.UTF_8);
		final String[] inputHeader=new String[15];Arrays.fill(inputHeader,"unused");inputHeader[1]="assembly";inputHeader[3]="taxid";inputHeader[6]="fasta";
		final Path inputs=write(dir,"inputs.tsv",String.join("\t",inputHeader)+"\n"+input("N",nuclear,hash)+"\n"+input("M",mito,"NA")+"\n");
		final Path verified=write(dir,"verified.tsv",hash+"\t"+nuclear+"\n"+hash+"\t"+mito+"\n");
		final Path models=write(dir,"models.fa",">mitochondria_insect\n"+repeat("ACGT",5)+"\n");
		final Path truth=write(dir,"truth.tsv",EukPanelTruthFreeze.HEADER+"\n"+truth("N","chr1",15,19,"+","BELOW_GA")+"\n"+truth("M","mito1",20,30,"+","ADMITTED")+"\n");
		final Path history=write(dir,"history.tsv","assembly\tnative_seqid\tstrand\tstart\tend\tfold_family\thistory_class\tmask_eligible\tsource_gff\tsource_sha80\tsource_line\toriginal_type\toriginal_attributes\tmask_reason\nN\tchr1\t-\t11\t12\tSSU\tHIST_ONLY\t1\tsource.gff\t"+hash+"\t1\trRNA\tproduct=test\tfixture\n");
		final String[] mh=new String[33],mr=new String[33];Arrays.fill(mh,"unused");Arrays.fill(mr,"NA");mh[31]="mask_eligible";mr[1]="N";mr[8]="chr1";mr[17]="15";mr[18]="16";mr[20]="+";mr[31]="0";
		final Path masks=write(dir,"masks.tsv","# fixture source masks\n"+String.join("\t",mh)+"\n"+String.join("\t",mr)+"\n");
		final ObjectMap<String,String> flags=new ObjectMap<String,String>(String.class,String.class);
		flags.put("verified",verified.toString());
		final String tree=dir.resolve("tree.tsv").toString();NcrnaMitoCrossLibraryTest.tree(tree);flags.put("tree",tree);
		flags.put("audit",auditPath.toString());flags.put("inputs",inputs.toString());flags.put("models",models.toString());flags.put("truth",truth.toString());flags.put("history",history.toString());flags.put("masks",masks.toString());
		flags.put("p2baseline",p2base.toString());flags.put("p2a",p2a.toString());flags.put("p2b",p2b.toString());flags.put("p3baseline",p3base.toString());flags.put("p3a",p3a.toString());flags.put("p3b",p3b.toString());
		final String prefix=dir.resolve("result").toString();flags.put("outprefix",prefix);
		new NcrnaMitoUnmatchedEvidence().run(flags,1,1);
		final java.util.List<String> rows=Files.readAllLines(Paths.get(prefix+".cases.tsv"),StandardCharsets.UTF_8);check(rows.size()==3,"Exactly one new unmatched call and one mitochondrial truth retained");
		final String[] n=rows.get(1).split("\t",-1),m=rows.get(2).split("\t",-1);
		check(n[1].equals("NUCLEAR_NEW_UNMATCHED") && n[4].equals("10") && n[8].equals("mitochondria_insect"),"Old unmatched coordinate must not become a new case");
		check(n[10].equals("NA") && n[11].equals("ORIGINAL_VERIFIER_NOT_RECORDED") && n[25].equals("UNRESOLVED"),"No identity or biological-origin inference from rounded scores or compartment labels");
		check(n[12].equals("80") && n[13].equals("80") && n[14].equals("40") && Double.parseDouble(n[15])==.5,"GC must use the whole original contig and explicit ACGT denominator");
		check(n[17].equals("2") && n[18].equals("1") && n[28].equals("chr1 original nuclear scaffold"),"BelowGA and eligible historical masks retained; original header preserved");
		check(m[19].equals("1") && m[20].equals("0") && m[21].equals("1") && m[22].equals("FN"),"All RNA overlaps, including opposite strand, and explicit zero arms survive without relabeling truth");
		check(m[31].equals("0") && m[32].equals("1") && m[33].equals("0") && m[35].equals("1"),"An arm with no RNA still retains separate CDS and other-feature overlaps at the organelle truth locus");
		check(n.length==42 && rows.get(0).split("\t",-1).length==42 && n[37].equals("4") && n[38].equals("2:Fixture phylum") && n[39].equals("3:Fixture class"),"Host ranks must follow the original roster taxid through the pinned tree, appended without shifting old columns");
		check(m[37].equals("999") && m[38].equals("MISSING_TAXID") && m[39].equals("MISSING_TAXID") && m[25].equals("UNRESOLVED"),"Missing taxonomy must remain explicit without dropping a locus or inferring mitochondrial origin");
		check(m[27].equals(hash) && m[40].equals("NA") && m[41].equals(verified.toString()) && n[40].equals(hash),"Resolve only original P3 NA digest from an exact verified source binding, preserving both provenances");
		check(NcrnaMitoUnmatchedEvidence.overlap(10,20,20,30)==1 && NcrnaMitoUnmatchedEvidence.overlap(10,20,21,30)==0,"Inclusive overlap distinguishes touching from adjacent intervals");
		flags.put("outprefix",dir.resolve("wrong_count").toString());boolean rejected=false;
		try{new NcrnaMitoUnmatchedEvidence().run(flags,2,1);}catch(IllegalArgumentException e){rejected=true;}
		check(rejected,"Unexpected selected populations must fail before writing a partial evidence table");
		final Path wrongBinding=write(dir,"wrong_verified.tsv",hash+"\t"+nuclear+"\n"+hash+"\t"+mito+".wrong\n");
		flags.put("verified",wrongBinding.toString());flags.put("outprefix",dir.resolve("wrong_binding").toString());rejected=false;
		try{new NcrnaMitoUnmatchedEvidence().run(flags,1,1);}catch(IllegalArgumentException e){rejected=e.getMessage().contains("exact FASTA path");}
		check(rejected,"Matching assembly alone must not resolve a digest from a different FASTA path");
		final Path conflictingBinding=write(dir,"conflicting_verified.tsv","00000000000000000002\t"+nuclear+"\n"+hash+"\t"+mito+"\n");
		flags.put("verified",conflictingBinding.toString());flags.put("outprefix",dir.resolve("conflicting_binding").toString());rejected=false;
		try{new NcrnaMitoUnmatchedEvidence().run(flags,1,1);}catch(IllegalArgumentException e){rejected=e.getMessage().contains("disagrees with original roster");}
		check(rejected,"An explicit ledger must never replace a conflicting non-NA roster digest");
		System.out.println("MITO_UNMATCHED_EVIDENCE_TEST_PASS original_rows whole_contig_GC masks opposite_strand zero_RNA_CDS_other no_origin_guess denominator");
	}
	private static String audit(String arm,String panel,String assembly,String contig,String compartment,String role,String id,String label,int start,int stop,String strand,String outcome){
		return String.join("\t",arm,panel,assembly,contig,compartment,"SSU",role,id,label,Integer.toString(start),Integer.toString(stop),strand,outcome,"NA","NA","NA","EDGE","FROZEN_LOCUS_TRUTH_V2");
	}
	private static String input(String id,Path fasta,String hash){final String[] row=new String[15];Arrays.fill(row,"NA");row[1]=id;row[3]=id.equals("N")?"4":"999";row[6]=fasta.toString();row[7]=hash;return String.join("\t",row);}
	private static String truth(String assembly,String contig,int start,int stop,String strand,String status){final String[] row=new String[35];Arrays.fill(row,"NA");row[0]="1";row[1]=assembly;row[2]=contig;row[3]=strand;row[4]=Integer.toString(start);row[5]=Integer.toString(stop);row[26]=status;row[32]="TRUTH_V2";return String.join("\t",row);}
	private static String gff(String contig,int start,int stop,String strand,String family,String model){return String.join("\t",contig,"fixture","rRNA",Integer.toString(start),Integer.toString(stop),"10.00",strand,"0","model:"+model+",family:"+family);}
	private static String repeat(String value,int n){final StringBuilder b=new StringBuilder();for(int i=0;i<n;i++){b.append(value);}return b.toString();}
	private static Path write(Path dir,String name,String text)throws Exception{final Path path=dir.resolve(name);Files.write(path,text.getBytes(StandardCharsets.UTF_8));return path;}
	private static void check(boolean ok,String why){if(!ok){throw new AssertionError(why);}}
}
