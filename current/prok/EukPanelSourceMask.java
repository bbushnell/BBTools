package prok;

import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import map.ObjectMap;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Preserve source-only, outside-target SRP annotations as provenance-bound masks.
 * This is not truth admission. The producer pins the input table and source manifest.
 * G11 2026-10-18: mask at SRP union, same strand, reciprocal overlap >=50% of both,
 * after admitted-locus matching; duplicate calls on consumed admitted loci remain FP.
 * @author Ganyu */
public final class EukPanelSourceMask {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args,EukPanelSourceMask.class,false);
		String in=null,sources=null,out=null;long expected=-1;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0,"Expected flag=value");
			final String k=arg.substring(0,eq).toLowerCase(Locale.ROOT),v=arg.substring(eq+1);
			if(k.equals("in")){in=v;}else if(k.equals("sources")){sources=v;}
			else if(k.equals("outprefix")){out=v;}else if(k.equals("expectedrows")){expected=Long.parseLong(v);}
			else{throw new IllegalArgumentException("Unknown flag: "+k);}
		}
		require(in!=null && sources!=null && out!=null && expected>=0,"Explicit input, source binding, output and row denominator required");
		Shared.setThreads(1);ReadWrite.setZipThreads(1);
		run(in,EukPanelGffValidation.hashes(sources),out,expected);Shared.closeStream(pp.outstream);
	}
	public static void run(String path,ObjectMap<String,String> hashes,String prefix,long expected){
		assert(hashes!=null && expected>=0):"Masks must refer to a pinned source manifest and exact review-table denominator";
		RrnaResourceIO.fresh(prefix,".out_of_map.tsv.gz",".mask_summary.tsv");
		final ByteFile in=RrnaResourceIO.open(path);
		final ByteStreamWriter out=RrnaResourceIO.writer(prefix+".out_of_map.tsv.gz"),summary=RrnaResourceIO.writer(prefix+".mask_summary.tsv");
		final LineParser1 p=new LineParser1('\t');final ByteBuilder b=new ByteBuilder();
		final ObjectMap<String,long[]> last=new ObjectMap<String,long[]>(String.class,long[].class);
		long rows=0,masks=0,eligible=0;
		try{
			RrnaResourceIO.header(in,HEADER,path);
			out.println("# OUT_OF_MAP mask: SRP biological union (small and large); same strand; reciprocal overlap >=0.5 of both intervals; apply after admitted 1:1 matching; duplicates of consumed admitted loci remain FP. No TP/FN admission. Only mask_eligible=1 rows may mask.");
			out.println(HEADER+"\tmask_id\tbiological_family_union\trfam_accessions\tmask_eligible\tmask_status");
			for(byte[] line;(line=in.nextLine())!=null;){
				p.set(line);require(p.terms()==28 && p.termEquals("UNRESOLVED_SOURCE_REVIEW_ONLY",27),"Expected original 28-column review evidence");rows++;
				final String assembly=p.parseString(1),source=p.parseString(2),hash=p.parseString(3);
				require(!assembly.isEmpty() && hash.equals(hashes.get(source)),"Source path/hash is absent from the pinned manifest: "+source);
				final String key=assembly+'\t'+source;long[] previous=last.get(key);
				if(previous==null){previous=new long[1];last.put(key,previous);}
				final long sourceLine=p.parseLong(4);require(sourceLine>previous[0],"Duplicate or unordered original source identity: "+key+":"+sourceLine);previous[0]=sourceLine;
				if(!p.termEquals("OUTSIDE_CURRENT_SUBTYPE_MAP",26)){continue;}
				final String attrs=p.parseString(22),type=p.parseString(16);
				require(type.endsWith("RNA") && !type.equals("mRNA") && EukPanelCmProposal.srpSubtype(attrs,true).equals("OUTSIDE_CURRENT_SUBTYPE_MAP"),"Outside-map mask must retain independently classified non-messenger SRP source evidence");
				final String accessions=outsideAccessions(attrs);
				final boolean valid=p.termEquals("VALID_INTERVAL",5);
				if(valid){require(p.termEquals("VALID_1_BASED_INCLUSIVE",7) && !p.termEquals("NA",8) && p.parseLong(17)>0 && p.parseLong(18)>=p.parseLong(17) && p.parseLong(18)<=p.parseLong(9),"Valid mask coordinates disagree with original contig bounds");}
				final boolean strand=p.termEquals("+",20) || p.termEquals("-",20),use=valid && strand;
				masks++;if(use){eligible++;}
				out.print(b.clear().append(line).tab().append(masks).tab().append("SRP").tab().append(accessions).tab().append(use?1:0)
					.tab().append(use?"OUT_OF_MAP":valid?"UNSCORABLE_STRAND":"UNSCORABLE_INTERVAL").nl());
			}
			require(rows==expected,"Source review row count changed: "+rows+" != "+expected);
			summary.println("review_rows\toutside_map_rows\teligible_masks\tunscorable_rows");
			summary.print(b.clear().append(rows).tab().append(masks).tab().append(eligible).tab().append(masks-eligible).nl());
		}finally{final boolean error=in.close();RrnaResourceIO.finish(out,summary);require(!error,"Source review read failed");}
		System.err.println("SOURCE_MASK_COMPLETE review_rows="+rows+" outside_map_rows="+masks+" eligible_masks="+eligible+" no_truth_admission=true");
	}
	static String outsideAccessions(String attrs){
		assert(attrs!=null):"Original attributes are retained; accession extraction must not guess subtype from length";
		final String db=EukPanelCmProposal.gffAttribute(attrs,"Dbxref");require(db!=null,"Outside-map source lacks RFAM provenance");
		final ByteBuilder b=new ByteBuilder();
		for(String token:db.split(",")){if(token.startsWith("RFAM:")){
			final String a=token.substring(5);require(a.matches("RF[0-9]{5}"),"Malformed RFAM accession: "+token);
			require(!a.equals("RF00169") && !a.equals("RF04183") && !a.equals("RF01854") && !a.equals("RF01857"),"A shipped SRP subtype cannot be an outside-map mask");
			if(b.length()>0){b.append(';');}b.append(a);
		}}
		require(b.length()>0,"Outside-map source requires at least one explicit RFAM accession");return b.toString();
	}
	public static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	public static final String HEADER=EukPanelCmProposal.SOURCE_HEADER+"\tsource_srp_subtype\tproposal_status";
}
