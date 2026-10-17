package prok;

import java.util.Arrays;
import java.util.HashSet;

import fileIO.ByteFile;
import parse.LineParser1;

/** Independent byte-preservation and accounting checks for a completed pilot.
 * @author Keqing
 */
public final class AnnotateGenesOutputTest {
	public static void main(String[] args){
		require(args.length==3, "Expected original GFF, annotated GFF and evidence TSV");
		final ByteFile before=ByteFile.makeByteFile(args[0], false), after=ByteFile.makeByteFile(args[1], false);
		final LineParser1 lp=new LineParser1('\t');
		long lines=0, cds=0, assigned=0, noHit=0, unsupported=0;
		boolean embedded=false;
		try{
			for(byte[] original=before.nextLine(); original!=null; original=before.nextLine()){
				lines++; final byte[] annotated=after.nextLine(); require(annotated!=null, "Annotated GFF lost line "+lines);
				if(Arrays.equals(original, "##FASTA".getBytes(java.nio.charset.StandardCharsets.US_ASCII))){embedded=true;}
				lp.set(original);
				if(embedded || original.length==0 || original[0]=='#' || lp.terms()!=9 || !lp.termEquals("CDS", 2)){
					require(Arrays.equals(original, annotated), "Changed non-CDS GFF line "+lines); continue;
				}
				cds++;
				final int prefix=original.length>=2 && original[original.length-1]=='.' && original[original.length-2]=='\t'
					? original.length-1 : original.length;
				require(annotated.length>prefix, "Missing appended CDS annotation on line "+lines);
				for(int i=0; i<prefix; i++){require(original[i]==annotated[i], "Changed existing CDS byte on line "+lines+", offset "+i);}
				lp.set(annotated); require(lp.terms()==9, "Annotation changed GFF column count on line "+lines);
				final String attributes=lp.parseString(8), status=CdsTranslator.attribute(attributes, "bbgap_status");
				if("assigned".equals(status)){
					assigned++; require(CdsTranslator.attribute(attributes, "bbgap_family_id")!=null, "Assigned CDS lacks family ID");
				}else if("no_hit".equals(status)){noHit++;}
				else if("unsupported".equals(status)){
					unsupported++; require(CdsTranslator.attribute(attributes, "bbgap_reason")!=null, "Unsupported CDS lacks reason");
				}else{throw new IllegalArgumentException("Missing or unknown annotation status on line "+lines);}
			}
			require(after.nextLine()==null, "Annotation added unexpected GFF rows");
		}finally{final boolean a=before.close(), b=after.close(); require(!a && !b, "GFF verification read failed");}
		final ByteFile evidence=ByteFile.makeByteFile(args[2], false);
		long rows=0, ea=0, en=0, eu=0, previous=0; boolean header=false; int width=0;
		final HashSet<String> packs=new HashSet<String>();
		try{
			for(byte[] line=evidence.nextLine(); line!=null; line=evidence.nextLine()){
				if(line.length==0){continue;} lp.set(line);
				if(line[0]=='#'){
					if(lp.termEquals("#pack", 0)){
						require(lp.terms()==6 && packs.add(lp.parseString(1)), "Malformed or duplicate pack evidence header");
					}
					continue;
				}
				if(!header){
					width=lp.terms();
					require((width==23 || (width==24 && lp.termEquals("pack_id", 23)))
						&& lp.termEquals("line", 0) && lp.termEquals("status", 9), "Unknown evidence header");
					require(width==23 ? packs.isEmpty() : !packs.isEmpty(), "Evidence pack headers and columns disagree");
					header=true; continue;
				}
				require(lp.terms()==width, "Unexpected evidence column count");
				rows++; final long index=lp.parseLong(0); require(index>previous && index<=lines, "Evidence rows are not in GFF order"); previous=index;
				if(lp.termEquals("assigned", 9)){ea++;}else if(lp.termEquals("no_hit", 9)){en++;}
				else if(lp.termEquals("unsupported", 9)){eu++;}else{throw new IllegalArgumentException("Unknown evidence status");}
				if(width==24){
					require(lp.termEquals("assigned", 9) ? packs.contains(lp.parseString(23)) : lp.termEquals(".", 23),
						"Assigned evidence must name a declared pack; non-hits must not name a winner");
				}
			}
		}finally{require(!evidence.close(), "Evidence read failed");}
		require(header && cds==rows && assigned==ea && noHit==en && unsupported==eu, "GFF/evidence CDS counters disagree");
		System.err.println("ANNOTATE_GENES_OUTPUT_PASS lines="+lines+" cds="+cds+" assigned="+assigned+" no_hit="+noHit+" unsupported="+unsupported+" original_bytes_preserved=true");
	}
	private static void require(boolean ok, String reason){if(!ok){throw new IllegalArgumentException(reason);}}
}
