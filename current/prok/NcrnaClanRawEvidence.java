package prok;

import java.util.ArrayList;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectMap;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Attributes every raw clan hit before deciding full-model support.
 * Weak additional hits and noncontainment geometry cannot hide a qualified hit.
 * Original-source rows and separate default/sensitive scan counts are retained.
 * @author Brian, Raiden
 */
public final class NcrnaClanRawEvidence {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args, NcrnaClanRawEvidence.class, false);
		String loci=null, models=null, out=null, family=null;int expected=-1;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0, "Use flag=value arguments");
			final String k=arg.substring(0, eq).toLowerCase(Locale.ROOT), v=arg.substring(eq+1);
			if(k.equals("loci")){loci=v;}else if(k.equals("rawmodels")){models=v;}else if(k.equals("outprefix")){out=v;}
			else if(k.equals("sourcefamily")){family=v;}else if(k.equals("expected")){expected=Integer.parseInt(v);}
			else{throw new IllegalArgumentException("Unknown option: "+k);}
		}
		require(loci!=null && models!=null && out!=null && family!=null && family.matches("RF[0-9]{5}") && expected>0,
			"Bind the original source ledger, all six raw scan tables and complete source denominator");
		Shared.setThreads(1);
		try{final NcrnaClanRawEvidence x=new NcrnaClanRawEvidence();x.load(loci, family, expected);x.models(models);x.write(out);}
		finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	private void load(String path, String family, int expected){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{
			RrnaResourceIO.header(in, NcrnaCmLocusAudit.HEADER, path);
			for(byte[] row; (row=in.nextLine())!=null; ){
				p.set(row);require(p.terms()==21, "Raw attribution needs the intact original21-field ledger");
				final String id=p.parseString(0);require(id.startsWith(family+"_") && !byId.contains(id), "Foreign or duplicated source ID: "+id);
				final Locus l=new Locus(id, row, p.parseInt(9), p.parseInt(18), p.parseInt(19));
				byId.put(id, l);loci.add(l);
			}
		}finally{require(!in.close(), "Original source read failed");}
		require(loci.size()==expected, "Original locus denominator changed");
	}
	private void models(String path){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');int rows=0;
		try{
			RrnaResourceIO.header(in, "model_id\tclen\tga_bits\tstage\thits", path);
			for(byte[] row; (row=in.nextLine())!=null; ){
				p.set(row);require(p.terms()==5, "Every raw table needs explicit model, criteria and scan stage");
				final String model=p.parseString(0), stage=p.parseString(3), hits=p.parseString(4);int m=-1;
				for(int i=0; i<IDS.length; i++){if(model.equals(IDS[i])){m=i;}}
				final int s=stage.equals("DEFAULT")?0:stage.equals("SENSITIVE")?1:-1;
				require(m>=0 && s>=0 && !seen[m][s] && p.parseInt(1)==CLEN[m] && p.parseDouble(2)==GA[m] && hits.startsWith("/"),
					"Pin one default and one sensitive table per Rfam15.1 clan model, without inherited or altered criteria");
				seen[m][s]=true;rows++;scan(m, s, hits);
			}
		}finally{require(!in.close(), "Raw model manifest read failed");}
		require(rows==6, "All three models need both independent scan stages");
	}
	private void scan(int m, int stage, String path){
		assert(m>=0 && m<3 && stage>=0 && stage<2):"Each hit belongs to one declared model and scan stage";
		final RrnaCmTruthBuilder.Criteria criteria=new RrnaCmTruthBuilder.Criteria(IDS[m], IDS[m], CLEN[m], CLEN[m], GA[m]);
		for(RrnaCmTruthBuilder.Hit h:RrnaCmTruthBuilder.readHits(new RrnaCmTruthBuilder.Genome("clan", "contexts", path, null), criteria)){
			final Locus l=byId.get(h.seqid);require(l!=null, "Raw hit target is absent from original source population: "+h.seqid);
			require(Math.max(h.seqFrom, h.seqTo)<=l.length+l.left+l.right, "Raw hit exceeds its original context: "+h.seqid);
			final Evidence e=l.evidence[m];e.raw[stage]++;
			if(h.strand!='+'){continue;}
			final int overlap=overlap(l, h);if(overlap==0){continue;}e.overlapping[stage]++;
			if(e.diagnostic==null || overlap>overlap(l, e.diagnostic)){e.diagnostic=h;}
			if(!h.status.equals("QUALIFIED")){continue;}
			e.qualified[stage]++;
			// A supported default-scan span remains the coordinate authority; sensitive
			// evidence fills missing support. Within a stage, maximize source overlap,
			// never compare bit scores between different covariance models.
			if(e.best==null || stage<e.stage || (stage==e.stage && overlap>overlap(l, e.best))){e.best=h;e.stage=stage;e.tied=false;}
			else if(stage==e.stage && overlap==overlap(l, e.best) && (h.seqFrom!=e.best.seqFrom || h.seqTo!=e.best.seqTo)){e.tied=true;}
		}
	}
	private void write(String out){
		RrnaResourceIO.fresh(out, ".models.tsv", ".selection.tsv", ".summary.tsv");
		for(String id:IDS){RrnaResourceIO.fresh(out, "."+id+".audit.tsv.gz");}
		// Refuse an ambiguous coordinate choice before opening output files.
		for(Locus l:loci){for(int m=0; m<3; m++){require(!l.evidence[m].tied, "Equal-overlap qualified cores require coordinate review: "+l.id+" model="+IDS[m]);}}
		final ByteStreamWriter manifest=RrnaResourceIO.writer(out+".models.tsv"), selection=RrnaResourceIO.writer(out+".selection.tsv"), summary=RrnaResourceIO.writer(out+".summary.tsv");
		final ByteBuilder b=new ByteBuilder();
		try{
			manifest.println("model_id\tclen\tga_bits\taudit");
			selection.println("locus_id\tmodel_id\tdefault_raw\tsensitive_raw\tdefault_overlapping\tsensitive_overlapping\tdefault_qualified\tsensitive_qualified\tselected_stage\tselected_table_line\tseq_from\tseq_to\tgeometry\tstatus");
			summary.println("model_id\tloci\tfull_supported\tgeometry_corrected\twith_multiple_qualified_events\tdefault_raw_events\tsensitive_raw_events");
			for(int m=0; m<3; m++){
				final String audit=out+'.'+IDS[m]+".audit.tsv.gz";final ByteStreamWriter rows=RrnaResourceIO.writer(audit);
				int supported=0, corrected=0, multiple=0, rawDefault=0, rawSensitive=0;
				try{
					rows.println(NcrnaCoreRecovery.AUDIT_HEADER);
					for(Locus l:loci){
						final Evidence e=l.evidence[m];final boolean full=e.best!=null;final RrnaCmTruthBuilder.Hit h=full?e.best:e.diagnostic;
						final int raw=e.raw[0]+e.raw[1], overlaps=e.overlapping[0]+e.overlapping[1];rawDefault+=e.raw[0];rawSensitive+=e.raw[1];
						if(full){supported++;if(h.seqFrom!=l.left+1 || h.seqTo!=l.left+l.length){corrected++;}if(e.qualified[0]+e.qualified[1]>1){multiple++;}}
						final String state=full?"CM_FULL_SPAN_ANY_HIT_SUPPORTED":h!=null?"RAW_NO_QUALIFYING_FULL_HIT":raw>0?"NO_OVERLAPPING_SENSE_HIT":"NO_REPORTED_HIT";
						b.clear().append(l.source).tab().append(state).tab().append(raw).tab().append(overlaps).tab();
						if(h==null){b.append("NA\tNA\tNA\tNA\tNA\tNA\t0\tNA\tNA\tNONE");}
						else{b.append(h.mdlFrom).tab().append(h.mdlTo).tab().append(h.seqFrom).tab().append(h.seqTo).tab().append(h.trunc)
							.tab().append(h.score).tab().append(overlap(l, h)).tab().append(h.seqFrom-l.left-1).tab().append(l.left+l.length-h.seqTo).tab().append(geometry(l, h));}
						rows.print(b.nl());
						b.clear().append(l.id).tab().append(IDS[m]).tab().append(e.raw[0]).tab().append(e.raw[1]).tab().append(e.overlapping[0]).tab().append(e.overlapping[1])
							.tab().append(e.qualified[0]).tab().append(e.qualified[1]).tab();
						if(full){b.append(e.stage==0?"DEFAULT":"SENSITIVE").tab().append(h.line).tab().append(h.seqFrom).tab().append(h.seqTo).tab().append(geometry(l, h)).append("\tFULL_SUPPORT");}
						else{b.append("NA\tNA\tNA\tNA\tNA\tUNRESOLVED");}selection.print(b.nl());
					}
				}finally{RrnaResourceIO.finish(rows);}
				manifest.print(b.clear().append(IDS[m]).tab().append(CLEN[m]).tab().append(GA[m], 2).tab().append(audit).nl());
				summary.print(b.clear().append(IDS[m]).tab().append(loci.size()).tab().append(supported).tab().append(corrected).tab().append(multiple).tab().append(rawDefault).tab().append(rawSensitive).nl());
			}
		}finally{RrnaResourceIO.finish(manifest, selection, summary);}
		System.err.println("NCRNA_CLAN_RAW_EVIDENCE_PASS loci="+loci.size()+" models=3 stages=2; raw event counts are not physical loci");
	}
	private static int overlap(Locus l, RrnaCmTruthBuilder.Hit h){return Math.max(0, Math.min(l.left+l.length, h.seqTo)-Math.max(l.left+1, h.seqFrom)+1);}
	private static String geometry(Locus l, RrnaCmTruthBuilder.Hit h){
		assert(overlap(l, h)>0):"Geometry only describes same-sense source-overlapping hits";
		final boolean hit=h.seqFrom>=l.left+1 && h.seqTo<=l.left+l.length, source=h.seqFrom<=l.left+1 && h.seqTo>=l.left+l.length;
		return hit?(source?"EXACT":"HIT_WITHIN_ANNOTATION"):source?"ANNOTATION_WITHIN_HIT":"PARTIAL_OVERLAP";
	}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------         Inner classes        ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Locus {
		Locus(String id_, byte[] source_, int length_, int left_, int right_){
			require(length_>0 && left_>=0 && right_>=0, "Original sense context must contain a nonempty source locus");
			id=id_;source=source_;length=length_;left=left_;right=right_;
			for(int i=0; i<3; i++){evidence[i]=new Evidence();}
		}
		final String id;final byte[] source;final int length, left, right;final Evidence[] evidence=new Evidence[3];
	}
	private static final class Evidence {
		final int[] raw=new int[2], overlapping=new int[2], qualified=new int[2];
		RrnaCmTruthBuilder.Hit best, diagnostic;int stage=-1;boolean tied;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Fields/constants       ----------------*/
	/*--------------------------------------------------------------*/

	private final ObjectMap<String,Locus> byId=new ObjectMap<String,Locus>(String.class, Locus.class);
	private final ArrayList<Locus> loci=new ArrayList<Locus>();private final boolean[][] seen=new boolean[3][2];
	private static final String[] IDS={"RF00174", "RF01482", "RF01689"};
	private static final int[] CLEN={189, 149, 81};private static final double[] GA={31.5, 50, 36.5};
}
