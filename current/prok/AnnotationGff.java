package prok;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashMap;

import fileIO.ByteFile;
import gff.GffLine;
import parse.LineParser1;
import structures.ByteBuilder;
import prot.ProteinFamilyAssigner;

/** Keeps original GFF bytes alongside the CDS translations used for annotation.
 * All output uses LF line endings. Existing attributes are never reserialized.
 * @author Keqing
 */
final class AnnotationGff {

	static AnnotationGff load(final String path, final HashMap<String,byte[]> scaffolds, final int table){
		final AnnotationGff doc=new AnnotationGff();
		final ByteFile input=ByteFile.makeByteFile(path, false);
		final LineParser1 fields=new LineParser1('\t');
		final HashMap<String,Cds> ids=new HashMap<String,Cds>(), parents=new HashMap<String,Cds>();
		int fallback=table;
		boolean embedded=false, features=false;
		try{
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				final int index=doc.lines.size(); doc.lines.add(line);
				if(embedded || line.length==0){continue;}
				if(line[0]=='#'){
					final String text=new String(line, StandardCharsets.UTF_8);
					if(text.equals("##FASTA")){embedded=true;}
					else if(text.startsWith("##transl_table ") || text.startsWith("##transl_table\t")){
						require(!features, "Translation-table directive follows features at GFF line "+(index+1));
						fallback=Integer.parseInt(text.substring(14).trim());
					}else if(text.startsWith("##genetic-code custom")){
						require(!features, "Custom-code directive follows features at GFF line "+(index+1)); fallback=-1;
					}
					continue;
				}
				features=true; fields.set(line);
				require(fields.terms()==9, "Expected nine GFF columns at line "+(index+1));
				for(int i=0; i<9; i++){require(!fields.termEquals("", i), "Empty GFF field at line "+(index+1)+", column "+(i+1));}
				final GffLine feature=new GffLine(line);
				final byte[] bases=scaffolds.get(decode(feature.seqid));
				require(bases!=null, "Missing FASTA contig '"+feature.seqid+"' at GFF line "+(index+1));
				require(feature.start>=1 && feature.stop>=feature.start && feature.stop<=bases.length,
					"GFF interval outside contig at line "+(index+1)+": "+feature.seqid+":"+feature.start+"-"+feature.stop);
				rejectOwned(feature.attributes, index+1);
				if(!"CDS".equals(feature.type)){continue;}
				final Cds cds=new Cds(index, feature, CdsTranslator.translate(feature, bases, fallback));
				doc.cds.add(cds);
				markRepeated(ids, CdsTranslator.attribute(feature.attributes, "ID"), cds);
				final String parent=CdsTranslator.attribute(feature.attributes, "Parent");
				if(parent!=null){
					//Multiple CDS segments for a transcript require splicing, outside C1.
					for(String token : parent.split(",", -1)){
						require(!token.isEmpty(), "Empty CDS Parent at GFF line "+(index+1));
						markRepeated(parents, token, cds);
					}
				}
			}
		}finally{require(!input.close(), "Failed reading GFF: "+path);}
		assert(doc.cds.size()<=doc.lines.size()) : "Each contiguous CDS must refer to one retained GFF row";
		return doc;
	}

	private static void markRepeated(final HashMap<String,Cds> seen, final String key, final Cds cds){
		if(key==null){return;}
		final Cds prior=seen.putIfAbsent(decode(key), cds);
		if(prior!=null && prior!=cds){prior.reason=cds.reason="compound_cds";}
	}

	/** Existing owned fields require an explicit replacement policy, not duplicates. */
	private static void rejectOwned(final String attributes, final int line){
		if(attributes==null){return;}
		for(int start=0; start<attributes.length(); ){
			int end=attributes.indexOf(';', start); if(end<0){end=attributes.length();}
			if(attributes.startsWith("bbgap_", start)){
				throw new IllegalArgumentException("Input already has BBGAP attributes at GFF line "+line);
			}
			start=end+1;
		}
	}

	/** Appends a first owned attribute, replacing only GFF's empty-attribute dot. */
	static ByteBuilder beginAnnotation(final ByteBuilder out, final byte[] raw){
		assert(raw!=null && raw.length>0) : "Annotation must extend an existing parsed GFF feature";
		out.clear();
		if(raw.length>=2 && raw[raw.length-1]=='.' && raw[raw.length-2]=='\t'){
			out.append(raw, 0, raw.length-1);
		}else{
			out.append(raw); if(raw[raw.length-1]!=';'){out.append(';');}
		}
		return out;
	}

	/** Percent-encodes UTF-8 and reserved separators in newly appended values only. */
	static ByteBuilder escaped(final ByteBuilder out, final String value){
		assert(value!=null) : "GFF attribute values must be explicit, never Java null";
		for(byte b : value.getBytes(StandardCharsets.UTF_8)){
			final int c=b&255;
			if(c<=32 || c>=127 || c=='%' || c==';' || c=='=' || c=='&' || c==','){
				out.append('%').append(HEX[c>>>4]).append(HEX[c&15]);
			}else{out.append(b);}
		}
		return out;
	}

	/** Decodes GFF identifiers for matching FASTA names; '+' is a literal plus. */
	static String decode(final String value){
		if(value.indexOf('%')<0){return value;}
		final byte[] raw=value.getBytes(StandardCharsets.UTF_8); final ByteBuilder out=new ByteBuilder(raw.length);
		for(int i=0; i<raw.length; i++){
			if(raw[i]!='%'){out.append(raw[i]); continue;}
			require(i+2<raw.length, "Truncated GFF escape: "+value);
			final int a=Character.digit((char)raw[++i], 16), b=Character.digit((char)raw[++i], 16);
			require(a>=0 && b>=0, "Invalid GFF escape: "+value); out.append((byte)((a<<4)|b));
		}
		return new String(out.toBytes(), StandardCharsets.UTF_8);
	}

	static final class Cds {
		Cds(int line_, GffLine feature_, CdsTranslator.Translation translation_){
			line=line_; feature=feature_; translation=translation_; reason=translation.reason;
		}
		String status(){return reason!=null ? "unsupported" : hit.assigned() ? "assigned" : "no_hit";}
		final int line;
		final GffLine feature;
		final CdsTranslator.Translation translation;
		String reason;
		ProteinFamilyAssigner.Hit hit;
	}
	private static void require(boolean ok, String reason){if(!ok){throw new IllegalArgumentException(reason);}}
	final ArrayList<byte[]> lines=new ArrayList<byte[]>();
	final ArrayList<Cds> cds=new ArrayList<Cds>();
	private static final char[] HEX="0123456789ABCDEF".toCharArray();
}
