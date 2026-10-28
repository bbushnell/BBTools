package prok;

import java.util.ArrayList;
import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import stream.Writer;
import stream.WriterFactory;
import structures.ByteBuilder;

/** Unions three DEFAULT CM evidence tables and preserves unresolved contexts
 * byte-for-byte for a sensitive followup. The upstream support flag requires
 * full-model/GA/no-truncation plus reciprocal-half call overlap; it is not a
 * biological verdict or evidence that excluded contexts lack nearby full hits.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaReviewCmUnion {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args, NcrnaReviewCmUnion.class, false);
		String loci=null, contexts=null, models=null, out=null;int expected=-1;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0, "Use flag=value");
			final String k=arg.substring(0, eq).toLowerCase(Locale.ROOT), v=arg.substring(eq+1);
			if(k.equals("loci")){loci=v;}else if(k.equals("in")){contexts=v;}else if(k.equals("models")){models=v;}
			else if(k.equals("outprefix")){out=v;}else if(k.equals("expected")){expected=Integer.parseInt(v);}else{throw new IllegalArgumentException(k);}
		}
		require(loci!=null && contexts!=null && models!=null && out!=null && expected>0, "Bind original review metadata/FASTA and three default evidence tables");
		Shared.setThreads(1);Shared.TRIM_READ_DESCRIPTION=false;
		try{final NcrnaReviewCmUnion x=new NcrnaReviewCmUnion();x.loci(loci, expected);x.models(models);x.write(contexts, out);}
		finally{Shared.closeStream(pp.outstream);}
	}
	private void loci(String path, int expected){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{
			RrnaResourceIO.header(in, NcrnaReviewUnion.LOCUS_HEADER, path);
			for(byte[] row; (row=in.nextLine())!=null; ){
				p.set(row);require(p.terms()==14 && p.termEquals("cobalamin", 5) && p.termEquals("PREDICTION_UNRESOLVED", 13), "Only the original B12 prediction union is eligible");
				final String id=p.parseString(0);final Row r=new Row(id, row.clone(), p.parseInt(10), p.parseInt(11), p.parseInt(12));
				require(!id.isEmpty() && !byId.contains(id) && r.left>=0 && r.right>=0 && r.length==p.parseInt(8)-p.parseInt(7)+1+r.left+r.right,
					"Review IDs and original clipped context geometry must remain exact");byId.put(id, r);rows.add(r);
			}
		}finally{require(!in.close(), "Review metadata read failed");}
		require(rows.size()==expected, "The complete original review population is required");
	}
	private void models(String path){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');final boolean[] seen=new boolean[3];
		try{
			RrnaResourceIO.header(in, "model_id\tstage\tloci", path);
			for(byte[] row; (row=in.nextLine())!=null; ){
				p.set(row);require(p.terms()==3 && p.termEquals("DEFAULT", 1), "Only default-pipeline qualification can select the subsequent sensitive subset");
				int index=-1;for(int i=0; i<IDS.length; i++){if(p.termEquals(IDS[i], 0)){index=i;}}
				require(index>=0 && !seen[index], "One evidence table per pinned B12 clan model is required");seen[index]=true;evidence(p.parseString(2), index);
			}
		}finally{require(!in.close(), "Model evidence manifest read failed");}
		for(boolean s:seen){require(s, "All three clan models must be considered before a context is unresolved");}
	}
	private void evidence(String path, int model){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');final ObjectSet<String> seen=new ObjectSet<String>(String.class);
		try{
			RrnaResourceIO.header(in, NcrnaReviewJoin.CM_HEADER, path);
			for(byte[] line; (line=in.nextLine())!=null; ){
				p.set(line);require(p.terms()==5, "Native per-locus CM evidence schema changed");
				final String id=p.parseString(0);final Row r=byId.get(id);final int raw=p.parseInt(1), overlap=p.parseInt(2), full=p.parseInt(3);
				require(r!=null && seen.add(id) && raw>=overlap && overlap>=full && full>=0, "Each model must account for each original context exactly once");
				require(p.termEquals(full==0?"NO_FULL_GA_OVERLAPPING_HIT":full==1?"ONE_FULL_GA_OVERLAPPING_HIT":"MULTIPLE_FULL_GA_OVERLAPPING_HITS", 4), "CM status must agree with native qualification count");
				r.raw+=raw;r.overlap+=overlap;if(full>0){r.mask|=1<<model;}
			}
		}finally{require(!in.close(), "Model locus evidence read failed");}
		require(seen.size()==rows.size(), "Missing model evidence cannot be treated as no support");
	}
	private void write(String fasta, String out){
		RrnaResourceIO.fresh(out, ".loci.tsv", ".summary.tsv", ".unresolved.loci.tsv", ".unresolved.contexts.fa.gz");
		final ObjectSet<String> seen=new ObjectSet<String>(String.class);final ArrayList<Read> selected=new ArrayList<Read>();
		RrnaResourceIO.stream(fasta, r->{
			final String id=RrnaResourceIO.first(r.id);final Row row=byId.get(id);
			require(row!=null && seen.add(id) && r.bases!=null && r.length()==row.length, "FASTA must preserve every original review ID and context length");
			require(r.id.equals(id+" lflank="+row.left+" rflank="+row.right), "Preserve the original sense-oriented flank header; never recompute it from a model hit");
			if(row.mask==0){selected.add(new Read(r.bases.clone(), null, r.id, selected.size(), false));}
		});
		require(seen.size()==rows.size(), "Every original context must be observed before writing a sensitive subset");
		final ByteStreamWriter all=RrnaResourceIO.writer(out+".loci.tsv"), subset=RrnaResourceIO.writer(out+".unresolved.loci.tsv"), summary=RrnaResourceIO.writer(out+".summary.tsv");
		final Writer sequence=WriterFactory.makeWriter(FileFormat.testOutput(out+".unresolved.contexts.fa.gz", FileFormat.FASTA, null, true, false, false, true), null, 1);
		sequence.start();final ByteBuilder b=new ByteBuilder();int supported=0, unresolved=0;long bases=0;
		try{
			all.println("review_id\tfull_model_mask_RF00174_RF01482_RF01689\traw_hit_events\tsense_overlap_events\tstatus");subset.println(NcrnaReviewUnion.LOCUS_HEADER);
			for(Row r:rows){
				if(r.mask==0){unresolved++;bases+=r.length;subset.println(r.rawMetadata);}else{supported++;}
				all.print(b.clear().append(r.id).tab().append(r.mask).tab().append(r.raw).tab().append(r.overlap).tab()
					.append(r.mask==0?"UNRESOLVED_FOR_SENSITIVE_SCAN":"HAS_DEFAULT_FULL_GA_RECIPROCAL_HALF_NOT_ADJUDICATED").nl());
			}
			require(supported+unresolved==rows.size() && unresolved==selected.size(), "Evidence partitions all contexts without merging alternative predictions");
			sequence.add(selected, 0);
			summary.println("original_contexts\twith_default_full_ga_reciprocal_half_support\tunresolved_contexts\tunresolved_context_bases");
			summary.print(b.clear().append(rows.size()).tab().append(supported).tab().append(unresolved).tab().append(bases).nl());
		}finally{final boolean failed=sequence.poisonAndWait();RrnaResourceIO.finish(all, subset, summary);require(!failed, "Sensitive subset FASTA must finish successfully");}
		System.err.println("NCRNA_REVIEW_CM_UNION_PASS contexts="+rows.size()+" supported="+supported+" unresolved="+unresolved+"; support is not a biological verdict");
	}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final class Row {
		Row(String i, byte[] bytes, int l, int r, int n){id=i;rawMetadata=bytes;left=l;right=r;length=n;}
		final String id;final byte[] rawMetadata;final int left, right, length;int mask;long raw, overlap;
	}
	private final ObjectMap<String,Row> byId=new ObjectMap<String,Row>(String.class, Row.class);
	private final ArrayList<Row> rows=new ArrayList<Row>();
	private static final String[] IDS={"RF00174", "RF01482", "RF01689"};
}
