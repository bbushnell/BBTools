package test.prok;

import prok.*;

import java.io.File;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;

/** Synthetic byte/schema fixtures; placeholder hashes are deliberately not integrity claims.
 * @author Ganyu */
public final class EukPanelCmEvidenceTest {
	public static void main(String[] args){
		require(args.length==1 && args[0].startsWith("outdir="),"Expected outdir=fresh_directory");
		final String root=args[0].substring(7);require(new File(root).mkdir(),"Fresh private fixture required");
		final String hash="00000000000000000000", fasta=root+"/source.fa", gff=root+"/source.gff", cm=root+"/models.cm";
		write(fasta,">ctg\nACGT\n");write(gff,"# Synthetic source pointer only\n");write(cm,"synthetic models; metadata is explicit\n");
		write(root+"/command.txt","synthetic cmsearch recipe\n");write(root+"/manifest.sha80","synthetic receipt\n");write(root+"/cmsearch","synthetic binary pointer\n");
		write(root+"/sources.tsv",hash+'\t'+fasta+'\n'+hash+'\t'+gff+'\n');
		write(root+"/inputs.tsv","compartment\tassembly\tbin\ttaxid\tbases\tsequence_records\tfasta\tselected_fasta_sha80\tgff\tannotation_status\tgbff\tgbff_record_accession\tmetadata_status\tannotation_source_fasta\tpairing_provenance\n"
			+"nuclear\tA\ttest\t1\t1000\t1\t"+fasta+"\tNA\t"+gff+"\tBOUND\tNA\tNA\tTEST\tNA\tTEST\n");
		write(root+"/contigs.tsv","selection_group\tassembly\tseqid\tlength\tcompartment\tevidence\tsource\tsource_sha80\tmatched_accession\traw_location\tassembly_unit\tsource_organelle_qualifier\tnative_taxid\tassembly_taxid\ttaxonomy_status\n"
			+"nuclear\tA\tctg\t1000\tNUCLEAR\tASSEMBLY_REPORT_REFSEQ\tsynthetic\t"+hash+"\tctg\tChromosome\tPrimary Assembly\tNA\t0\t1\tNOT_TAGGED\n");
		final String models=EukPanelCmEvidence.MODEL_HEADER+'\n'
			+model(hash,"RF00001","5S",120,38,cm)+model(hash,"RF00002","5.8S",154,42,cm)+model(hash,"RF00005","tRNA",71,29,cm);
		write(root+"/models.tsv",models);
		final String run="r1\tA\t"+root+"/raw.tbl\t"+hash+'\t'+cm+'\t'+hash+'\t'+fasta+'\t'+hash+'\t'+root+"/cmsearch\t"+hash+'\t'+root+"/command.txt\t"+hash+'\t'+root+"/manifest.sha80\t"+fasta+'\t'+hash+"\tINFERNAL_1_1_5_TBLOUT17\n";
		write(root+"/runs.tsv",EukPanelCmEvidence.RUN_HEADER+'\n'+run);
		final String raw="# Program: cmsearch\n# Version: 1.1.5 (synthetic)\n"
			+hit("ctg","RF00001",10,30,"+",42,"!")+hit("ctg","RF00002",20,25,"+",40,"?")
			+hit("ctg","RF00001",60,40,"-",45,"!")+hit("orphan","RF00001",5,10,"+",41,"!")
			+hit("ctg","RF00001",999,1005,"+",39,"!")+"# [ok]\n";
		write(root+"/raw.tbl",raw);
		write(root+"/features.tsv",EukPanelCmEvidence.FEATURE_HEADER+'\n'
			+feature(gff,hash,2,"VALID_INTERVAL","10","30","+")
			+feature(gff,hash,3,"VALID_INTERVAL","40","60","-")
			+feature(gff,hash,4,"VALID_INTERVAL","200","210","+")
			+feature(gff,hash,5,"OUT_OF_RANGE","995","1020","+")
			+feature(gff,hash,6,"INVALID_COORDINATES","oops","30","+")
			+feature(gff,hash,7,"VALID_INTERVAL","10","30","."));
		final String[] files={"inputs.tsv","sources.tsv","contigs.tsv","features.tsv","models.tsv","runs.tsv","raw.tbl","models.cm","source.fa","source.gff","command.txt","manifest.sha80","cmsearch","missing_footer.tbl","bad_runs.tsv"};
		final StringBuilder verified=new StringBuilder();for(String f:files){verified.append(hash).append('\t').append(root).append('/').append(f).append('\n');}
		write(root+"/verified.tsv",verified.toString());
		run(root,"runs.tsv","pass");verify(root);
		write(root+"/missing_footer.tbl",raw.replace("# [ok]\n",""));
		write(root+"/bad_runs.tsv",EukPanelCmEvidence.RUN_HEADER+'\n'+run.replace("/raw.tbl","/missing_footer.tbl"));
		expectFailure(root,"bad_runs.tsv","missing_footer","Incomplete");
		write(root+"/duplicate_runs.tsv",EukPanelCmEvidence.RUN_HEADER+'\n'+run+run);
		// A separate metadata manifest is needed even for a deliberate bad fixture.
		write(root+"/duplicate_verified.tsv",verified.toString()+hash+'\t'+root+"/duplicate_runs.tsv\n");
		boolean rejected=false;try{run(root,"duplicate_runs.tsv","duplicate","duplicate_verified.tsv");}catch(IllegalArgumentException e){require(e.getMessage().contains("Duplicate run"),"Wrong duplicate failure: "+e);rejected=true;}
		require(rejected,"Duplicate selected scan must fail rather than double evidence");
		write(root+"/extra_models.tsv",models+model(hash,"RF00003","other",100,30,cm));
		write(root+"/extra_verified.tsv",verified.toString()+hash+'\t'+root+"/extra_models.tsv\n");
		rejected=false;try{run(root,"runs.tsv","extra_roster","extra_verified.tsv","extra_models.tsv");}
		catch(IllegalArgumentException e){require(e.getMessage().contains("Exact global accession roster"),"Wrong model-roster failure: "+e);rejected=true;}
		require(rejected,"An extra accession must fail rather than substitute a same-size model set");
		System.out.println("PASS all raw hits, no-hit model, below-GA inclusion distinction, cross-model overlap, minus strand, unresolved/out-of-range hits, original-source conservation and anomaly retention; incomplete/duplicate scans rejected");
	}
	static String model(String hash,String acc,String family,int clen,int ga,String path){return hash+'\t'+acc+'\t'+family+'\t'+acc+'\t'+clen+'\t'+ga+'\t'+path+'\t'+hash+'\n';}
	static String hit(String target,String model,int from,int to,String strand,int score,String inc){return target+" - name "+model+" cm 1 "+(Math.abs(to-from)+1)+' '+from+' '+to+' '+strand+" no 1 0.5 0.0 "+score+" 1e-12 "+inc+" description\n";}
	static String feature(String gff,String hash,int line,String status,String from,String to,String strand){return "nuclear\tA\t"+gff+'\t'+hash+'\t'+line+'\t'+status+"\tEXACT_LITERAL\t"+status+"\tctg\t1000\tNUCLEAR\t0\t0\tNOT_TAGGED\tctg\tsynthetic\tncRNA\t"+from+'\t'+to+"\t.\t"+strand+"\t.\tID=f"+line+'\n';}
	static void run(String root,String runs,String out){run(root,runs,out,"verified.tsv");}
	static void run(String root,String runs,String out,String verified){
		run(root,runs,out,verified,"models.tsv");
	}
	static void run(String root,String runs,String out,String verified,String models){
		EukPanelCmEvidence.main(new String[]{"inputs="+root+"/inputs.tsv","sources="+root+"/sources.tsv","contigs="+root+"/contigs.tsv","contigsources="+root+"/sources.tsv","features="+root+"/features.tsv","models="+root+'/'+models,"runs="+root+'/'+runs,"verified="+root+'/'+verified,"outprefix="+root+'/'+out,"expected=1","expectedmodels=3"});
	}
	static void verify(String root){
		final String[] expected={"MULTIPLE_MODELS","ONE_HIT","NO_SAME_STRAND_OVERLAP","NOT_ASSOCIATED_OUT_OF_RANGE","NOT_ASSOCIATED_INVALID_COORDINATES","STRAND_UNSPECIFIED"};
		final ByteFile sources=RrnaResourceIO.open(root+"/pass.sources.tsv.gz"), original=RrnaResourceIO.open(root+"/features.tsv");final LineParser1 p=new LineParser1('\t');int n=0;
		try{sources.nextLine();original.nextLine();for(byte[] line;(line=sources.nextLine())!=null;){
			p.set(line);require(p.terms()==26 && n<expected.length && p.termEquals(expected[n],23),"Source status/count/order mismatch");
			final byte[] old=original.nextLine();require(old!=null && line.length>old.length && line[old.length]=='\t',"Original source row missing");
			for(int i=0;i<old.length;i++){require(old[i]==line[i],"All original source fields must be unchanged");}n++;
		}require(n==6 && original.nextLine()==null,"Conserve every source row");}finally{require(!sources.close() && !original.close(),"Source fixture read failed");}
		final ByteFile summary=RrnaResourceIO.open(root+"/pass.summary.tsv");n=0;
		try{summary.nextLine();for(byte[] line;(line=summary.nextLine())!=null;){p.set(line);require(p.terms()==13,"Summary schema");
			if(p.termEquals("RF00001",2)){require(p.parseLong(4)==4 && p.parseLong(5)==2 && p.parseLong(6)==1 && p.parseLong(7)==1 && p.parseLong(12)==2,"5S raw/mapped/anomaly/overlap conservation");}
			else if(p.termEquals("RF00002",2)){require(p.parseLong(4)==1 && p.parseLong(8)==0 && p.parseLong(11)==0 && p.parseLong(12)==1,"Retain below-GA ? hit without accepting it");}
			else{require(p.termEquals("RF00005",2) && p.parseLong(4)==0,"No-hit model must retain its denominator");}n++;
		}require(n==3,"All model scan summaries retained");}finally{require(!summary.close(),"Summary fixture read failed");}
		final ByteFile raw=RrnaResourceIO.open(root+"/pass.hits.tsv.gz");n=0;
		final String[] ranges={"VALID_INTERVAL","VALID_INTERVAL","VALID_INTERVAL","UNRESOLVED_SEQID","OUT_OF_RANGE"};
		try{p.set(raw.nextLine());require(p.terms()==29,"Raw hit header width");for(byte[] line;(line=raw.nextLine())!=null;){
			p.set(line);require(p.terms()==29 && n<5 && p.parseLong(4)==n+3 && p.termEquals(ranges[n],13),"Every raw hit retains its original line and explicit interval status");
			if(n==1){require(p.termEquals("RF00002",6) && p.termEquals("?",23) && p.parseDouble(27)==-2,"Below-GA raw hit must remain unchanged");}
			if(n==2){require(p.termEquals("-",21) && p.parseInt(19)==60 && p.parseInt(20)==40,"Retain original minus-oriented sequence coordinates");}
			n++;
		}require(n==5,"Five raw rows must survive");}finally{require(!raw.close(),"Raw hit fixture read failed");}
		final ByteFile edges=RrnaResourceIO.open(root+"/pass.associations.tsv.gz");n=0;
		try{p.set(edges.nextLine());require(p.terms()==10,"Association header width");for(byte[] line;(line=edges.nextLine())!=null;){
			p.set(line);require(p.terms()==10,"Association row width");
			if(n<3){require(p.termEquals("OVERLAP",3) && p.termEquals("r1",4),"Expected real physical overlap");
				final long expectedSource=n<2?2:3, expectedHit=n+3, expectedOverlap=n==1?6:21;
				require(p.parseLong(2)==expectedSource && p.parseLong(5)==expectedHit && p.parseLong(7)==expectedOverlap,"Exact plus/minus source/hit overlap identity");
				require(Math.abs(p.parseDouble(8)-(n==1?6.0/21:1))<1e-8 && p.parseDouble(9)==1,"Reciprocal overlap fractions must use inclusive spans");
			}else{require(p.termEquals(expected[n-1],3) && p.termEquals("NA",7),"Zero-overlap/anomaly sentinel must retain reason rather than fabricate a numeric overlap");}
			n++;
		}require(n==7,"Three overlaps plus four source-status sentinels");}finally{require(!edges.close(),"Association fixture read failed");}
		final ByteFile cross=RrnaResourceIO.open(root+"/pass.hit_overlaps.tsv.gz");n=0;
		try{p.set(cross.nextLine());require(p.terms()==11,"Cross-hit header has eleven declared fields");for(byte[] line;(line=cross.nextLine())!=null;){
			p.set(line);require(p.terms()==11 && p.termEquals("A",0) && p.termEquals("ctg",1) && p.termEquals("+",2)
				&& p.parseLong(4)==3 && p.termEquals("RF00001",5) && p.parseLong(7)==4 && p.termEquals("RF00002",8)
				&& p.parseLong(9)==6 && p.termEquals("CROSS_MODEL",10),"Exact cross-model pair, strand and physical overlap required");n++;
		}require(n==1,"Exactly one cross-model pair; minus/unresolved/out-of-range hits are not fabricated overlaps");}finally{require(!cross.close(),"Cross-model fixture read failed");}
	}
	static void expectFailure(String root,String runs,String out,String reason){boolean failed=false;try{run(root,runs,out);}catch(IllegalArgumentException e){require(e.getMessage().contains(reason),"Wrong failure: "+e);failed=true;}require(failed,"Expected "+reason);}
	static void write(String path,String text){final ByteStreamWriter out=RrnaResourceIO.writer(path);out.print(text);RrnaResourceIO.finish(out);}
	static void require(boolean ok,String why){EukPanelCmEvidence.require(ok,why);}
}
