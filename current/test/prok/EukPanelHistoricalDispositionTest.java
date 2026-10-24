package test.prok;

import java.io.File;
import java.util.Arrays;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import prok.EukPanelHistoricalDisposition;
import prok.EukPanelHistoricalMask;
import prok.EukPanelTruthDisposition;
import prok.EukPanelTruthFreeze;
import prok.RrnaResourceIO;
import structures.ByteBuilder;

/** End-to-end historical decisions, strand-qualified baseline joins and conservation.
 * @author Brian, Raiden */
public final class EukPanelHistoricalDispositionTest {
	public static void main(String[] args){
		require(args.length==1 && args[0].startsWith("outdir="), "Fresh fixture directory required");
		final String dir=args[0].substring(7);require(new File(dir).mkdirs(), "Fixture directory must be fresh");
		final String truth=dir+"/truth.tsv", models=dir+"/models.tsv", evidence=dir+"/evidence.tsv", history=dir+"/history.tsv";
		write(truth, EukPanelTruthFreeze.HEADER, new String[][]{locus("L1", 500, "+", true), locus("L2", 300, "+", false), locus("L3", 700, "-", true)});
		write(models, "bundle_sha80\taccession\tfamily\tsubtype\tclen\tga\tsource_cm\tsource_cm_sha80",
			new String[][]{{HASH,"RF00005","tRNA","tRNA","100","50","fixture.cm",HASH}});
		final String[][] rows={row("H1",100,"SUPPORTED_GA",false),row("H2",100,"SUPPORTED_GA",false),
			row("H3",105,"SUPPORTED_GA",false),row("H4",300,"SUPPORTED_SUBGA",false),row("H5",500,"SUPPORTED_GA",false),
			row("H6",700,"UNSUPPORTED",true),row("H7",900,"UNSUPPORTED",true),row("H8",1100,"UNSUPPORTED",false)};
		write(evidence, EukPanelHistoricalDisposition.INPUT_HEADER, rows);
		final String[][] masks=new String[rows.length][];for(int i=0; i<rows.length; i++){masks[i]=Arrays.copyOfRange(rows[i], 1, 15);}
		write(history, EukPanelHistoricalMask.OUTPUT_HEADER, masks);
		EukPanelHistoricalDisposition.run(truth,evidence,models,history,HASH,dir+"/ok",3,8,8);
		final ByteFile in=RrnaResourceIO.open(dir+"/ok.dispositions.tsv.gz");final LineParser1 p=new LineParser1('\t');int n=0;
		try{
			require(in.nextLine()!=null, "Output header missing");
			for(byte[] line;(line=in.nextLine())!=null;){p.set(line);n++;require(p.terms()==122,"Disposition schema changed");
				if(n<=3){require(p.termEquals("ADMITTED",113) && p.termEquals("NONE",117) && p.termEquals("NONE",118),"New GA evidence must not acquire a baseline link");}
				if(n<=2){require(p.termEquals("HIST_V2_H1",114) && p.parseInt(115)==2 && p.termEquals("HIST_V2_H3",116),"Exact duplicates group, shifted overlapping evidence stays explicit");}
				if(n==4){require(p.termEquals("BELOW_GA",113) && p.termEquals("L2",118),"Below-GA existing locus must be visible, not duplicated or admitted");}
				if(n==5){require(p.termEquals("L1",117),"GA evidence must link existing admitted locus");}
				if(n==6){require(p.termEquals("ORIENTATION_ERROR_NO_NEW_LOCUS",113) && p.termEquals("L3",119),"Opposite evidence must join using the HIT strand");}
				if(n==7){require(p.termEquals("ORIENTATION_ERROR_NO_NEW_LOCUS",113) && p.termEquals("NONE",119),"Opposite support alone must not invent an existing counterpart");}
				if(n==8){require(p.termEquals("LIFT_MASK",113),"Unsupported evidence lifts its old mask");}
			}
		}finally{require(!in.close(),"Fixture output read failed");}require(n==8,"Every annotation must survive");
		final ByteFile delta=RrnaResourceIO.open(dir+"/ok.delta.tsv");try{delta.nextLine();p.set(delta.nextLine());
			final long[] want={8,4,1,3,2,1,1,4,2,3};for(int i=0; i<want.length; i++){require(p.parseLong(i+1)==want[i],"Wrong family delta field "+i);}
			require(delta.nextLine()==null,"Unexpected extra family");
		}finally{require(!delta.close(),"Delta read failed");}
		negative(truth,models,history,rows,dir,0,23,"SUPPORTED_SUBGA","class");
		negative(truth,models,history,rows,dir,0,29,"ffffffffffffffffffff","hash");
		negative(truth,models,history,rows,dir,0,83,"-","strand");
		negative(truth,models,history,rows,dir,0,13,"changed attributes","history");
		negative(truth,models,history,rows,dir,0,86,"99","model_length");
		final String[][] bad=rows.clone();bad[1]=rows[0];final String duplicate=dir+"/duplicate.tsv";
		write(duplicate,EukPanelHistoricalDisposition.INPUT_HEADER,bad);
		expectFailure(()->EukPanelHistoricalDisposition.run(truth,duplicate,models,history,HASH,dir+"/bad_duplicate",3,8,8));
		expectFailure(()->EukPanelHistoricalDisposition.run(truth,evidence,models,history,HASH,dir+"/bad_count",3,8,7));
		materialization(truth,models,rows,dir);
		System.out.println("HISTORICAL_DISPOSITION_TEST_PASS classes=true baseline_below_ga=true opposite_strand=true duplicates=true overlap_flags=true negatives=7");
	}
	static void materialization(String truth,String models,String[][] rows,String dir){
		assert(rows.length==8):"The materialization fixture reuses the already-checked evidence cases";
		final String[][] selected={rows[0],rows[3],rows[4],rows[5],rows[6],rows[7],row("H9",1300,"SUPPORTED_SUBGA",false)};
		selected[6][1]="empty";
		final String evidence=dir+"/materialize_evidence.tsv",history=dir+"/materialize_history.tsv";
		final String inputs=dir+"/panel_inputs.tsv";
		write(inputs,"compartment\tassembly\tbin\ttaxid\tbases\tsequence_records\tfasta\tselected_fasta_sha80\tgff\tannotation_status\tgbff\tgbff_record_accession\tmetadata_status\tannotation_source_fasta\tpairing_provenance",
			new String[][]{{"nuclear","assembly","0","1","2000","1","source.fa",HASH,"source.gff","fixture","NA","NA","fixture","source.fa","fixture"},
			{"nuclear","empty","nuclear_other","2","2000","1","source.fa",HASH,"source.gff","fixture","NA","NA","fixture","source.fa","fixture"}});
		write(evidence,EukPanelHistoricalDisposition.INPUT_HEADER,selected);
		final String[][] masks=new String[selected.length][];for(int i=0; i<selected.length; i++){masks[i]=Arrays.copyOfRange(selected[i],1,15);}
		write(history,EukPanelHistoricalMask.OUTPUT_HEADER,masks);
		EukPanelHistoricalDisposition.run(truth,evidence,models,history,HASH,dir+"/materialized",3,7,7,true,inputs,2);
		final ByteFile original=RrnaResourceIO.open(truth),combined=RrnaResourceIO.open(dir+"/materialized.truth_candidate.tsv.gz");
		final LineParser1 p=new LineParser1('\t'),q=new LineParser1('\t');int n=0;
		try{
			original.nextLine();combined.nextLine();
			for(byte[] line;(line=combined.nextLine())!=null;){p.set(line);n++;require(p.terms()==35,"Combined truth schema must remain35 columns");
				require(p.termEquals(EukPanelTruthFreeze.V2_CANDIDATE,32) && p.termEquals(EukPanelTruthFreeze.V2_POLICY,34),"Materialization cannot freeze itself");
				if(n<=3){q.set(original.nextLine());for(int i=0;i<35;i++){if(i!=32 && i!=34){require(p.parseString(i).equals(q.parseString(i)),"Baseline biological field changed at "+i);}}}
				else{require(p.termEquals(n==4?"HIST_V2_H1":"HIST_V2_H9",0) && p.termEquals(n==4?"ADMITTED":"BELOW_GA",26),"Only new GA/subGA loci may be appended");
					require(p.parseInt(18)==1 && p.termEquals("FULLSPAN_INTACT",27) && p.parseInt(31)==(n==4?1:0),"New locus geometry/member/boundary eligibility incorrect");
					if(n==5){require(p.termEquals("empty",1) && p.termEquals("nuclear_other",7),"An assembly absent from baseline must retain its roster bin");}}
			}require(n==5 && original.nextLine()==null,"Existing admitted/below-GA links must not duplicate baseline loci");
		}finally{require(!original.close() && !combined.close(),"Materialized truth read failed");}
		final ByteFile inventory=RrnaResourceIO.open(dir+"/materialized.historical_masks.tsv.gz");n=0;int retained=0;
		try{inventory.nextLine();for(byte[] line;(line=inventory.nextLine())!=null;){p.set(line);n++;require(p.terms()==14,"History schema changed");
			if(p.termEquals("1",7)){retained++;require(p.termEquals("HIST_LINKED_UNSCORED",6),"Only subGA historical masks remain");}
			if(p.parseInt(10)==6 || p.parseInt(10)==7){require(p.termEquals("ORIENTATION_ERROR",6) && p.termEquals("0",7),"Opposite-only support must lift original-strand mask without adding truth");}
		}}finally{require(!inventory.close(),"Materialized history read failed");}require(n==7 && retained==2,"Historical masks must conserve rows and retain only subGA");
		final ByteFile totals=RrnaResourceIO.open(dir+"/materialized.materialization.tsv");try{totals.nextLine();p.set(totals.nextLine());
			final long[] want={3,1,1,2,2,5,5};for(int i=0;i<want.length;i++){require(p.parseLong(i+1)==want[i],"Wrong materialized count at "+i);}
			require(totals.nextLine()==null,"Unexpected materialization family");
		}finally{require(!totals.close(),"Materialization totals read failed");}
		expectFailure(()->EukPanelHistoricalDisposition.run(truth,dir+"/evidence.tsv",models,dir+"/history.tsv",HASH,dir+"/ambiguous_materialize",3,8,8,true,inputs,2));
		final String[][] upgrade=selected.clone();upgrade[1]=row("H4",300,"SUPPORTED_GA",false);
		write(dir+"/upgrade_evidence.tsv",EukPanelHistoricalDisposition.INPUT_HEADER,upgrade);
		expectFailure(()->EukPanelHistoricalDisposition.run(truth,dir+"/upgrade_evidence.tsv",models,history,HASH,dir+"/upgrade_materialize",3,7,7,true,inputs,2));
		System.out.println("HISTORICAL_MATERIALIZATION_TEST_PASS baseline_preserved=true linked_not_duplicated=true masks=true ambiguous_and_promotion_rejected=true");
	}
	static void negative(String truth,String models,String history,String[][] original,String dir,int row,int field,String value,String name){
		assert(row>=0 && row<original.length):"Negative fixture changes exactly one specified input field";
		final String[][] changed=original.clone();changed[row]=original[row].clone();changed[row][field]=value;
		final String path=dir+"/negative_"+name+".tsv";write(path,EukPanelHistoricalDisposition.INPUT_HEADER,changed);
		expectFailure(()->EukPanelHistoricalDisposition.run(truth,path,models,history,HASH,dir+"/bad_"+name,3,8,8));
	}
	static void expectFailure(Runnable task){boolean caught=false;try{task.run();}catch(IllegalArgumentException e){caught=true;}require(caught,"Malformed input was silently accepted");}
	static String[] locus(String id,int start,String strand,boolean ga){
		final String[] a=new String[35];Arrays.fill(a,"0");a[0]=id;a[1]="assembly";a[2]="contig";a[3]=strand;a[4]=""+start;a[5]=""+(start+99);
		a[6]="tRNA";a[14]=ga?"1":"0";a[26]=ga?"ADMITTED":"BELOW_GA";a[32]="TRUTH_V1";a[34]=EukPanelTruthDisposition.LEGACY_POLICY;return a;
	}
	static String[] row(String id,int start,String classification,boolean opposite){
		final String[] a=new String[113];Arrays.fill(a,"NA");
		a[0]=id;a[1]="assembly";a[2]="contig";a[3]="+";a[4]=""+start;a[5]=""+(start+99);a[6]="tRNA";
		a[7]="HIST_ONLY";a[8]="1";a[9]="source.gff";a[10]=HASH;a[11]=id.substring(1);a[12]="tRNA";a[13]="product=tRNA";a[14]="NO_RAW_HIT_LINK";
		a[15]="2000";a[16]="NUCLEAR";a[17]="1";a[18]="2000";a[19]="source.fa";a[20]="100";a[21]="tRNA";a[22]=id;
		a[23]=classification;a[24]="PROPOSE_ADDITION_REQUIRES_ADJUDICATION";a[25]="fixture";a[26]=opposite?"1":"0";a[27]="NONE";a[28]=id+"/hits.tbl";a[29]=HASH;a[30]="TRUTH_V2_CANDIDATE_ADDENDUM_NOT_FROZEN";
		if(!classification.equals("UNSUPPORTED")){hit(a,79,start,"+",classification.equals("SUPPORTED_GA")?60:49);}
		if(opposite){hit(a,95,start,"-",60);}a[111]="1";a[112]=opposite?"1":"0";return a;
	}
	static void hit(String[] a,int c,int start,String strand,int bits){
		final String[] h={"RF00005","1",""+start,""+(start+99),strand,"1","100","100","no",""+bits,bits>=50?"2.6e-08":"1e-3","50",""+(bits-50),"100","1","1"};
		System.arraycopy(h,0,a,c,h.length);
	}
	static void write(String path,String header,String[][] rows){
		assert(!new File(path).exists()):"Fixture inputs must not overwrite previous evidence";
		final ByteStreamWriter w=RrnaResourceIO.writer(path);final ByteBuilder b=new ByteBuilder();
		try{w.println(header);for(String[] row:rows){b.clear();for(int i=0; i<row.length; i++){if(i>0){b.tab();}b.append(row[i]);}w.print(b.nl());}}
		finally{RrnaResourceIO.finish(w);}
	}
	static void require(boolean ok,String why){if(!ok){throw new AssertionError(why);}}
	static final String HASH="0123456789abcdef0123";
}
