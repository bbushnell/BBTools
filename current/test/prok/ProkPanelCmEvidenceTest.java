package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import prok.EukPanelCmEvidence;
import prok.RrnaResourceIO;

/** Exercises explicit prok mode on actual native-validator output, including duplicate seqids across assemblies.
 * @author Brian Bushnell, Raiden */
public final class ProkPanelCmEvidenceTest {
	public static void main(String[] args){
		check(args.length==1 && args[0].startsWith("outdir="),"Use outdir=fresh_directory");
		final String d=args[0].substring(7);check(new File(d).mkdir(),"Fresh fixture directory required");
		NcrnaFreshPanelValidationTest.main(new String[]{"outdir="+d+"/native"});
		final String n=d+"/native", cm=d+"/models.cm", executable=d+"/cmsearch", command=d+"/command.txt", seal=d+"/seal.sha80";
		write(cm,"Synthetic model pointer; metadata binds the two declared accessions.\n");write(executable,"Synthetic executable pointer.\n");write(command,"Synthetic pinned recipe.\n");write(seal,"Synthetic scan seal pointer.\n");
		write(d+"/models.tsv",EukPanelCmEvidence.MODEL_HEADER+'\n'
			+EukPanelCmEvidenceTest.model(HASH,"RF00001","5S",120,38,cm)
			+EukPanelCmEvidenceTest.model(HASH,"RF00005","tRNA",71,29,cm));
		final String header="# Program: cmsearch\n# Version: 1.1.5 (synthetic)\n";
		write(d+"/a.tbl",header+EukPanelCmEvidenceTest.hit("shared","RF00001",1,8,"+",42,"!")
			+EukPanelCmEvidenceTest.hit("tagged","RF00005",4,1,"-",32,"!")+"# [ok]\n");
		write(d+"/b.tbl",header+EukPanelCmEvidenceTest.hit("shared","RF00001",1,5,"+",42,"!")+"# [ok]\n");
		write(d+"/runs.tsv",EukPanelCmEvidence.RUN_HEADER+'\n'
			+runRow("a","GCF_1.1",d+"/a.tbl",cm,n+"/a.fa",executable,command,seal)
			+runRow("b","GCF_2.1",d+"/b.tbl",cm,n+"/b.fa",executable,command,seal));
		final String[] paths={n+"/good.assemblies.tsv",n+"/good.contigs.tsv",n+"/good.features.tsv.gz",n+"/sources.tsv",
			n+"/a.fa",n+"/b.fa",n+"/a.gff",n+"/b.gff",cm,executable,command,seal,d+"/a.tbl",d+"/b.tbl",d+"/runs.tsv",d+"/models.tsv",
			d+"/bad.contigs.tsv",d+"/bad.inputs.tsv",d+"/bad.features.tsv"};
		final StringBuilder verified=new StringBuilder();for(String path:paths){verified.append(HASH).append('\t').append(path).append('\n');}write(d+"/verified.tsv",verified.toString());
		run(d,"pass",n+"/good.assemblies.tsv",n+"/good.contigs.tsv",n+"/good.features.tsv.gz",true);
		final ArrayList<byte[]> hits=read(d+"/pass.hits.tsv.gz");check(hits.size()==4,"Every raw hit remains present");
		final LineParser1 p=new LineParser1('\t');
		p.set(hits.get(1));check(p.termEquals("GCF_1.1",0) && p.termEquals("PROKARYOTIC",11) && p.termEquals("VALID_INTERVAL",13) && p.parseLong(28)==2,
			"The first assembly's shared contig supports both preserved duplicate source rows");
		p.set(hits.get(2));check(p.termEquals("tid|11|tagged",10) && p.termEquals("-",21) && p.parseLong(28)==1,
			"Native-tag alias and minus-strand association survive prok adaptation");
		p.set(hits.get(3));check(p.termEquals("GCF_2.1",0) && p.termEquals("OUT_OF_RANGE",13) && p.parseLong(28)==0,
			"An identically named contig in another assembly must use its own four-base limit");
		final ArrayList<byte[]> source=read(d+"/pass.sources.tsv.gz"), original=read(n+"/good.features.tsv.gz");
		check(source.size()==original.size() && source.size()==8,"No source row or anomaly may disappear");
		for(int row=1;row<source.size();row++){
			final byte[] a=source.get(row), b=original.get(row);check(a.length>b.length && a[b.length]=='\t',"Original feature prefix missing");
			for(int i=0;i<b.length;i++){check(a[i]==b[i],"Prok reducer changed an original source field");}
		}
		p.set(source.get(7));check(p.termEquals("archaea",0) && p.termEquals("NOT_ASSOCIATED_OUT_OF_RANGE",23),"Archaeal identity and source anomaly must remain explicit");
		final ArrayList<byte[]> summary=read(d+"/pass.summary.tsv");check(summary.size()==5,"Both models need a row for both assemblies, including no hits");
		p.set(summary.get(4));check(p.termEquals("GCF_2.1",0) && p.termEquals("RF00005",2) && p.parseLong(4)==0,"No-hit tRNA denominator retained");
		final String contigs=text(n+"/good.contigs.tsv"), inputs=text(n+"/good.assemblies.tsv"), features=text(n+"/good.features.tsv.gz");
		write(d+"/bad.contigs.tsv",contigs.replace("PROKARYOTIC","NUCLEAR"));bad(d,"nuclear",n+"/good.assemblies.tsv",d+"/bad.contigs.tsv",n+"/good.features.tsv.gz",true,"not nuclear labels");
		write(d+"/bad.contigs.tsv",contigs.replace(n+"/a.fa\t"+HASH,n+"/b.fa\t"+HASH));bad(d,"cross_source",n+"/good.assemblies.tsv",d+"/bad.contigs.tsv",n+"/good.features.tsv.gz",true,"own assembly FASTA/hash");
		write(d+"/bad.contigs.tsv",contigs.replace("shared\t8\t","shared\t7\t"));bad(d,"length",n+"/good.assemblies.tsv",d+"/bad.contigs.tsv",n+"/good.features.tsv.gz",true,"conserve every selected assembly");
		write(d+"/bad.inputs.tsv",inputs.replace("\t12\t12\t2\t","\t13\t12\t2\t"));bad(d,"selected_bases",d+"/bad.inputs.tsv",n+"/good.contigs.tsv",n+"/good.features.tsv.gz",true,"positive observed dimensions");
		write(d+"/bad.features.tsv",features.replace("\tPROKARYOTIC\t","\tNUCLEAR\t"));bad(d,"feature_compartment",n+"/good.assemblies.tsv",n+"/good.contigs.tsv",d+"/bad.features.tsv",true,"contradicts observed contig");
		bad(d,"explicit_mode_required",n+"/good.assemblies.tsv",n+"/good.contigs.tsv",n+"/good.features.tsv.gz",false,"corrected input binding schema");
		System.out.println("PROK_CM_EVIDENCE_TEST_PASS native_validator_schema assembly_scoped_coordinates original_rows_preserved zero_hit_models explicit_domain source_binding_rejection euk_default_unchanged");
	}
	private static String runRow(String id,String assembly,String table,String cm,String fasta,String exe,String command,String seal){
		return id+'\t'+assembly+'\t'+table+'\t'+HASH+'\t'+cm+'\t'+HASH+'\t'+fasta+'\t'+HASH+'\t'+exe+'\t'+HASH+'\t'+command+'\t'+HASH+'\t'+seal+'\t'+fasta+'\t'+HASH+"\tINFERNAL_1_1_5_TBLOUT17\n";
	}
	private static void run(String d,String label,String inputs,String contigs,String features,boolean prok){
		EukPanelCmEvidence.main(new String[]{"inputs="+inputs,"sources="+d+"/native/sources.tsv","contigs="+contigs,"contigsources="+d+"/native/sources.tsv",
			"features="+features,"runs="+d+"/runs.tsv","models="+d+"/models.tsv","verified="+d+"/verified.tsv","outprefix="+d+'/'+label,"expected=2","expectedmodels=2","prokinputs="+prok});
	}
	private static void bad(String d,String label,String inputs,String contigs,String features,boolean prok,String reason){
		boolean rejected=false;try{run(d,label,inputs,contigs,features,prok);}catch(IllegalArgumentException e){check(e.getMessage().contains(reason),"Wrong failure for "+label+": "+e);rejected=true;}
		check(rejected,"Bad metadata accepted: "+label);
	}
	private static ArrayList<byte[]> read(String path){
		final ArrayList<byte[]> rows=new ArrayList<byte[]>();final ByteFile in=RrnaResourceIO.open(path);
		try{for(byte[] line;(line=in.nextLine())!=null;){rows.add(line);}}finally{check(!in.close(),"Fixture read failed");}return rows;
	}
	private static String text(String path){final StringBuilder b=new StringBuilder();for(byte[] line:read(path)){b.append(new String(line,StandardCharsets.UTF_8)).append('\n');}return b.toString();}
	private static void write(String path,String text){final ByteStreamWriter out=new ByteStreamWriter(path,true,false,false);out.start();out.print(text);RrnaResourceIO.finish(out);}
	private static void check(boolean ok,String why){if(!ok){throw new AssertionError(why);}}
	private static final String HASH="00000000000000000000";
}
