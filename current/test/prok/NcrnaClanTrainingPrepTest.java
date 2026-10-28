package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import prok.NcrnaClanEvidence;
import prok.NcrnaClanOverlay;
import prok.NcrnaClanRawEvidence;
import prok.NcrnaClanTrainingPrep;
import prok.NcrnaCmLocusAudit;
import prok.RrnaResourceIO;
import structures.ByteBuilder;

/** Checks complete routing, biological labels, exact sequence conservation and gain attribution.
 * @author Brian, Raiden
 */
public final class NcrnaClanTrainingPrepTest {

	public static void main(String[] args){
		require(args.length==1, "Supply a fresh fixture directory");final File dir=new File(args[0]);
		require(!dir.exists() && dir.mkdirs(), "Fixture must preserve earlier evidence");final String b=dir.getAbsolutePath(), a=b+"/rf14";
		NcrnaClanRawEvidenceTest.main(new String[]{a});
		previous(a+"/overlay.overlay.tsv.gz", a+"/previous.tsv", "RF01482_1_0");
		final String s=b+"/rf16";require(new File(s).mkdir(), "Create a distinct sibling source fixture");
		write(s+"/source.tsv", NcrnaCmLocusAudit.HEADER, "RF01689_0_0\tg16\tGCF_16.1\t16\tbacteria\tP\tchr16\t+\t1001-1080\t80\t1\tf\tfalse\tfalse\tANNOTATION_FULL_COPY\t2000\tfalse\t0\t200\t200\t1");
		final ByteBuilder seq=new ByteBuilder();for(int i=0; i<480; i++){seq.append("AGTCCGTA".charAt(i%8));}
		write(s+"/contexts.fa", ">RF01689_0_0 source=g16 contig=chr16 strand=+ start=1001 stop=1080 lflank=200 rflank=200", seq.toString());
		final String[] ids={"RF00174", "RF01482", "RF01689"};final int[] clens={189,149,81};final double[] ga={31.5,50,36.5};
		final ByteStreamWriter manifest=RrnaResourceIO.writer(s+"/rawmodels.tsv");
		try{
			manifest.println("model_id\tclen\tga_bits\tstage\thits");
			for(int i=0; i<3; i++){for(String stage:new String[]{"DEFAULT", "SENSITIVE"}){
				final String path=s+'/'+ids[i]+'.'+stage+".tbl";
				final ByteStreamWriter w=RrnaResourceIO.writer(path);
				try{w.println("# Program: cmsearch");w.println("# Version: 1.1.5 (fixture)");
					if(i==2 && stage.equals("DEFAULT")){w.println("RF01689_0_0 - model RF01689 cm 1 81 198 270 + no 1 0.5 0.0 60 1e-10 ! fixture");}w.println("# [ok]");}
				finally{RrnaResourceIO.finish(w);}
				manifest.println(ids[i]+"\t"+clens[i]+"\t"+ga[i]+"\t"+stage+'\t'+path);
			}}
		}finally{RrnaResourceIO.finish(manifest);}
		NcrnaClanRawEvidence.main(new String[]{"loci="+s+"/source.tsv", "rawmodels="+s+"/rawmodels.tsv", "sourcefamily=RF01689", "expected=1", "outprefix="+s+"/raw"});
		NcrnaClanEvidence.main(new String[]{"models="+s+"/raw.models.tsv", "sourcefamily=RF01689", "expected=1", "modelcount=3", "outprefix="+s+"/join"});
		NcrnaClanOverlay.main(new String[]{"joined="+s+"/join.loci.tsv.gz", "contexts="+s+"/contexts.fa", "sourcefamily=RF01689", "expected=1", "minlen=66", "maxlen=134", "evidencestage=RAW_ADJUDICATED", "outprefix="+s+"/overlay"});
		previous(s+"/overlay.overlay.tsv.gz", s+"/previous.tsv", "NONE");
		write(b+"/base.fa", ">RF00174_base source=gbase family=RF00174 status=ANNOTATION_FULL_COPY", seq.toString().substring(0, 150));
		write(b+"/inputs.tsv", "family\tfasta\toverlay\tprevious_overlay\tjoined\tselection\trawmodels\texpected",
			"RF00174\t"+b+"/base.fa\tNA\tNA\tNA\tNA\tNA\t1",
			input("RF01482", a, 3), input("RF01689", s, 1));
		NcrnaClanTrainingPrep.main(new String[]{"inputs="+b+"/inputs.tsv", "outprefix="+b+"/train", "expected=5"});
		final ArrayList<String> summary=lines(b+"/train.summary.tsv"), members=lines(b+"/train.members.tsv"), gains=lines(b+"/train.gains.tsv");
		require(summary.size()==4 && summary.get(1).equals("RF00174\t1\t0\t0") && summary.get(2).equals("RF01482\t2\t0\t2")
			&& summary.get(3).equals("RF01689\t2\t1\t1"), "All five members route once; gain counts retain source family rather than operational destination");
		boolean cross=false, defaultGain=false, sensitiveGain=false;
		for(String row:members){if(row.startsWith("RF01482_1_0\t")){cross=row.contains("\tRF01482\tB12clade\tRF01689\tHIGHEST_BITS_OPERATIONAL_ONLY\t60.0\t120\t");}}
		for(String row:gains){if(row.startsWith("RF01482_0_0\t")){defaultGain=row.contains("\tRF01482\tRF01482\tDEFAULT\t") && row.contains("\tPARTIAL_OVERLAP\t");}
			if(row.startsWith("RF01482_2_0\t")){sensitiveGain=row.contains("\tRF01482\tRF01482\tSENSITIVE\t");}}
		require(members.size()==6 && gains.size()==4 && cross && defaultGain && sensitiveGain,
			"Operational higher-bits routing cannot relabel a clan member; native gain rows distinguish default geometry from sensitive support");
		final ArrayList<String> clan=lines(b+"/train.clan.fa"), rf1=lines(b+"/train.RF00174.fa"), rf14=lines(b+"/train.RF01482.fa"), rf16=lines(b+"/train.RF01689.fa");
		require(clan.size()==10 && rf1.size()==2 && rf14.size()==4 && rf16.size()==4, "Each grouped FASTA and complete clan FASTA must conserve all actual records");
		require(rf1.get(1).equals(lines(b+"/base.fa").get(1)) && rf14.get(1).equals(lines(a+"/overlay.admitted.fa").get(1))
			&& rf16.get(1).equals(lines(a+"/overlay.admitted.fa").get(3)) && rf16.get(3).equals(lines(s+"/overlay.admitted.fa").get(1)),
			"Routing must preserve exact accepted sequences, including corrected cross-labeled and sibling cores");
		boolean failed=false;try{NcrnaClanTrainingPrep.main(new String[]{"inputs="+b+"/inputs.tsv", "outprefix="+b+"/wrong", "expected=4"});}catch(IllegalArgumentException e){failed=e.getMessage().contains("combined training population");}
		require(failed && !new File(b+"/wrong.summary.tsv").exists(), "Wrong all-member denominator must fail before publishing training resources");
		final ArrayList<String> tied=lines(a+"/join.loci.tsv.gz");
		for(int i=1; i<tied.size(); i++){if(tied.get(i).startsWith("RF01482_1_0\t")){
			final String[] fields=tied.get(i).split("\t", -1);fields[34]="60";tied.set(i, String.join("\t", fields));}}
		write(a+"/routing_tie.tsv", tied.toArray(new String[tied.size()]));
		final ArrayList<String> tieInputs=lines(b+"/inputs.tsv");
		for(int i=0; i<tieInputs.size(); i++){tieInputs.set(i, tieInputs.get(i).replace(a+"/join.loci.tsv.gz", a+"/routing_tie.tsv"));}
		write(b+"/tie_inputs.tsv", tieInputs.toArray(new String[tieInputs.size()]));failed=false;
		try{NcrnaClanTrainingPrep.main(new String[]{"inputs="+b+"/tie_inputs.tsv", "outprefix="+b+"/routing_tie", "expected=5"});}
		catch(IllegalArgumentException e){failed=e.getMessage().contains("Equal-score cross-labeled routing");}
		require(failed && !new File(b+"/routing_tie.summary.tsv").exists(), "Equal-score operational routing must not invent a directory assignment");
		final ArrayList<String> changed=lines(a+"/previous.tsv");
		for(int i=1; i<changed.size(); i++){changed.set(i, changed.get(i).replace("\tg0\t", "\twrong\t"));}
		write(a+"/changed_previous.tsv", changed.toArray(new String[changed.size()]));
		final ArrayList<String> badInputs=lines(b+"/inputs.tsv");
		for(int i=0; i<badInputs.size(); i++){badInputs.set(i, badInputs.get(i).replace(a+"/previous.tsv", a+"/changed_previous.tsv"));}
		write(b+"/bad_inputs.tsv", badInputs.toArray(new String[badInputs.size()]));failed=false;
		try{NcrnaClanTrainingPrep.main(new String[]{"inputs="+b+"/bad_inputs.tsv", "outprefix="+b+"/changed", "expected=5"});}
		catch(IllegalArgumentException e){failed=e.getMessage().contains("Previous/current overlays");}
		require(failed && !new File(b+"/changed.summary.tsv").exists(), "Gain attribution cannot silently change original source identity");
		System.out.println("NCRNA_CLAN_TRAINING_PREP_TEST_PASS");
	}
	private static String input(String family, String dir, int expected){
		return family+'\t'+dir+"/overlay.admitted.fa\t"+dir+"/overlay.overlay.tsv.gz\t"+dir+"/previous.tsv\t"+dir+"/join.loci.tsv.gz\t"+dir+"/raw.selection.tsv\t"+dir+"/rawmodels.tsv\t"+expected;
	}
	private static void previous(String current, String out, String keep){
		final ArrayList<String> rows=lines(current);final String header=rows.get(0);final ByteStreamWriter w=RrnaResourceIO.writer(out);
		try{
			w.println(header.substring(0, header.lastIndexOf('\t')));
			for(int i=1; i<rows.size(); i++){
				final String[] f=rows.get(i).split("\t", -1);final StringBuilder row=new StringBuilder();
				for(int c=0; c<21; c++){if(c>0){row.append('\t');}row.append(f[c]);}
				if(f[0].equals(keep)){for(int c=21; c<32; c++){row.append('\t').append(f[c]);}}
				else{row.append('\t').append(f[21]).append("\tUNRESOLVED\tfalse\tNONE\tNA\tFIRST_PASS\tNA\tNA\tNA\tNA\tUNRESOLVED");}
				w.println(row.toString());
			}
		}finally{RrnaResourceIO.finish(w);}
	}
	private static void write(String path, String... rows){
		assert(rows.length>0):"Fixture files have explicit complete content";
		final ByteStreamWriter w=RrnaResourceIO.writer(path);try{for(String row:rows){w.println(row);}}finally{RrnaResourceIO.finish(w);}
	}
	private static ArrayList<String> lines(String path){
		final ByteFile in=RrnaResourceIO.open(path);final ArrayList<String> rows=new ArrayList<String>();
		try{for(byte[] row; (row=in.nextLine())!=null; ){rows.add(new String(row, StandardCharsets.UTF_8));}}
		finally{require(!in.close(), "Fixture output read failed");}return rows;
	}
	private static void require(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
}
