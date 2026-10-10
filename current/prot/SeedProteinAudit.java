package prot;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * Read-only audit of seed proteins against the protein encoder's residue contract.
 * Separates edge-only stop markers, internal stops, and unsupported residues;
 * audit mode never edits or excludes sequences. Explicit trimedges mode writes
 * a derived copy that repairs edge-only stop records and fails on other defects.
 * Counts and raw lengths
 * must match every pinned family-manifest entry before the final PASS is written.
 * @author Brian Bushnell, Keqing
 */
public final class SeedProteinAudit {
	public static void main(String[] args){
		try{
			Shared.setThreads(1); Shared.AMINO_IN=true;
			Read.VALIDATE_IN_CONSTRUCTOR=false;//Audit raw symbols without Read's case/junk rewriting.
			PreParser pp=new PreParser(args, SeedProteinAudit.class, false);
			try{run(HmmComparisonData.options(pp.args));}finally{Shared.closeStream(pp.outstream);}
		}catch(Throwable failure){failure.printStackTrace(); System.exit(1);}
	}

	private static void run(HashMap<String,String> o) throws Exception{
		for(String key : o.keySet()){
			if(!Arrays.asList("mode", "manifest", "manifestsha80", "source", "out").contains(key)){throw new IllegalArgumentException("Unknown audit parameter: "+key);}
		}
		final String mode=o.getOrDefault("mode", "audit");
		require(mode.equals("audit") || mode.equals("trimedges"), "Expected mode=audit or mode=trimedges");
		final boolean trim=mode.equals("trimedges");
		final String manifest=HmmComparisonData.required(o, "manifest"), pin=HmmComparisonData.required(o, "manifestsha80");
		DigestSuffix.requireSuffix(pin, "manifest pin"); require(DigestSuffix.file(manifest).equals(pin), "Manifest hash differs");
		final Path source=Paths.get(HmmComparisonData.required(o, "source")), out=Paths.get(HmmComparisonData.required(o, "out"));
		require(!Files.exists(out), "Audit output must be a fresh directory");
		final List<String[]> rows=HmmComparisonData.rows(manifest);
		require(!rows.isEmpty() && String.join("\t", rows.get(0)).equals("active_index\tfamily_id\trep_id\tassigned_count\tsequence_bytes\tfasta_file\tfasta_sha80"), "Unsupported audit manifest");
		Files.createDirectory(out);
		if(trim){Files.createDirectory(out.resolve("members"));}
		final ByteStreamWriter normalizedManifest=trim ? HmmComparisonData.writer(out.resolve("members/family_manifest.tsv").toString()) : null;
		if(trim){normalizedManifest.println(String.join("\t", rows.get(0)));}
		final ByteStreamWriter families=HmmComparisonData.writer(out.resolve("families.tsv").toString());
		final ByteStreamWriter rejected=HmmComparisonData.writer(out.resolve("affected_records.tsv").toString());
		families.println("active_index\trep_id\trecords\traw_residues\tvalid\tedge_stop_only\tinternal_stop\tunsupported\tinvalid_representative");
		rejected.println("active_index\tprotein_id\tcategory\tlength\tleading_stops\ttrailing_stops\tinternal_stops\tunsupported_residues");
		long total=0, rawTotal=0, normalizedTotal=0, valid=0, edge=0, internal=0, unsupported=0, badRepresentatives=0;
		final ByteBuilder fastaRow=new ByteBuilder();
		for(int rowIndex=1; rowIndex<rows.size(); rowIndex++){
			final String[] row=rows.get(rowIndex);
			require(row.length==7, "Audit manifest must have seven fields");
			final String input=source.resolve(row[5]).toString();
			require(row[5].equals("rank_"+Integer.parseInt(row[0])+".faa"), "Unexpected seed filename: "+row[5]);
			DigestSuffix.requireSuffix(row[6], "family input pin"); require(DigestSuffix.file(input).equals(row[6]), "Family input changed: "+row[0]);
			final Path derived=out.resolve("members").resolve(row[5]);
			final ByteStreamWriter normalized=trim ? HmmComparisonData.writer(derived.toString()) : null;
			final Streamer st=StreamerFactory.makeStreamer(FileFormat.testInput(input, FileFormat.FASTA, null, true, false), null, true, -1);
			long n=0, raw=0, normalizedRaw=0, good=0, edgeOnly=0, hasInternal=0, badResidue=0;
			boolean representativeSeen=false, invalidRepresentative=false;
			st.start();
			try{
				for(ListNum<Read> list=st.nextList(); list!=null && list.size()>0; list=st.nextList()){
					for(Read read : list){
						require(read.mate==null && read.bases!=null, "Expected an unpaired protein record");
						final byte[] b=read.bases;
						n++; raw+=b.length;
						int leading=0, trailing=0, stops=0, other=0;
						while(leading<b.length && b[leading]=='*'){leading++;}
						while(trailing<b.length-leading && b[b.length-1-trailing]=='*'){trailing++;}
						for(int i=leading; i<b.length-trailing; i++){
							if(b[i]=='*'){stops++;}
							else if(!legal(b[i]&255)){other++;}
						}
						final boolean normal=b.length>trailing && leading==0 && trailing<=1 && stops==0 && other==0;
						if(read.id.equals(row[2])){require(!representativeSeen, "Duplicate representative in audit family"); representativeSeen=true; invalidRepresentative=!normal;}
						int from=0, to=b.length;
						if(normal){good++;}
						else{
							final String category;
							if(other>0 || leading+trailing==b.length){badResidue++; category="UNSUPPORTED_OR_EMPTY";}
							else if(stops>0){hasInternal++; category="INTERNAL_STOP";}
							else{edgeOnly++; category="EDGE_STOP_ONLY";}
							rejected.print(row[0]).tab().print(read.id).tab().print(category).tab().print(b.length).tab()
								.print(leading).tab().print(trailing).tab().print(stops).tab().print(other).nl();
							if(trim){
								require(category.equals("EDGE_STOP_ONLY"), "Edge repair refuses "+category+" in "+read.id);
								from=leading; to=b.length-trailing;
							}
						}
						normalizedRaw+=to-from;
						if(trim){normalized.print(fastaRow.clear().append('>').append(read.id).nl().append(b, from, to-from).nl());}
					}
				}
			}finally{
				st.close();
				if(normalized!=null){HmmComparisonData.close(normalized);}
				require(!st.errorState(), "Audit FASTA input I/O error");
			}
			require(n==Long.parseLong(row[3]) && raw==Long.parseLong(row[4]) && representativeSeen,
				"Audit count/representative mismatch at active index "+row[0]);
			require(good+edgeOnly+hasInternal+badResidue==n, "Audit categories must partition every input record");
			if(trim){
				normalizedManifest.print(row[0]).tab().print(row[1]).tab().print(row[2]).tab().print(n).tab().print(normalizedRaw)
					.tab().print(row[5]).tab().print(DigestSuffix.file(derived.toString())).nl();
			}
			families.print(row[0]).tab().print(row[2]).tab().print(n).tab().print(raw).tab().print(good).tab()
				.print(edgeOnly).tab().print(hasInternal).tab().print(badResidue).tab().print(invalidRepresentative ? 1 : 0).nl();
			total+=n; rawTotal+=raw; normalizedTotal+=normalizedRaw; valid+=good; edge+=edgeOnly; internal+=hasInternal; unsupported+=badResidue;
			if(invalidRepresentative){badRepresentatives++;}
			if(rowIndex%100==0){System.err.println("SEED_AUDIT_PROGRESS families="+rowIndex+" records="+total+" affected="+(edge+internal+unsupported));}
		}
		HmmComparisonData.close(families); HmmComparisonData.close(rejected);
		if(normalizedManifest!=null){HmmComparisonData.close(normalizedManifest);}
		require(DigestSuffix.file(manifest).equals(pin), "Audit manifest changed during scan");
		ByteStreamWriter bw=HmmComparisonData.writer(out.resolve("summary.tsv").toString());
		bw.println("families\trecords\traw_residues\tvalid\tedge_stop_only\tinternal_stop\tunsupported\tinvalid_representatives");
		bw.print(rows.size()-1).tab().print(total).tab().print(rawTotal).tab().print(valid).tab().print(edge).tab().print(internal).tab().print(unsupported).tab().print(badRepresentatives).nl();
		HmmComparisonData.close(bw);
		if(trim){
			bw=HmmComparisonData.writer(out.resolve("normalization.tsv").toString());
			bw.println("records_preserved\trecords_repaired\tsource_raw_residues\tnormalized_raw_residues\tremoved_edge_markers");
			bw.print(total).tab().print(edge).tab().print(rawTotal).tab().print(normalizedTotal).tab().print(rawTotal-normalizedTotal).nl();
			HmmComparisonData.close(bw);
		}
		bw=HmmComparisonData.writer(out.resolve("PASS").toString()); bw.println("SEED_PROTEIN_AUDIT_PASS"); HmmComparisonData.close(bw);
		System.err.println("SEED_PROTEIN_AUDIT_PASS records="+total+" affected="+(edge+internal+unsupported));
	}

	/** Mirrors Blosum62.encodeResidue: twenty ordinary residues or X/B/Z/J, never U/O or stops. */
	private static boolean legal(int c){
		if(c<0 || c>=128){return false;}
		final byte n=AminoAcid.acidToNumber[c];
		if(n>=0 && n<=19){return true;}
		c=Character.toUpperCase(c);
		return c=='X' || c=='B' || c=='Z' || c=='J';
	}
	private static void require(boolean ok, String message){if(!ok){throw new IllegalArgumentException(message);}}
}
