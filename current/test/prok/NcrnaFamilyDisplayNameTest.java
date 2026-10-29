package test.prok;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import gff.GffLine;
import prok.NcrnaFamilyRegistry;
import prok.Orf;
import prok.ProkObject;

/** Tests scientific names at the real GFF emission and tab-delimited parsing boundary.
 * @author Brian, Raiden */
public final class NcrnaFamilyDisplayNameTest {
	public static void main(String[] args){
		final String[][] names={{"prok5S", "5S"}, {"rnasep", "RNase P"}, {"tmrna", "tmRNA"},
			{"fmn", "FMN"}, {"sixs_rf00013", "6S"}, {"sixs_rf01685", "6S-Flavo"},
			{"srp_small", "SRP RNA (small)"}, {"srp_large", "SRP RNA (large)"},
			{"16S", "16S"}, {"23S", "23S"}, {"sam1", "SAM-I"},
			{"euk5S", "euk5S"}, {"cobalamin_rf00174", "cobalamin_rf00174"}};
		final boolean oldParse=GffLine.parseAttributes;GffLine.parseAttributes=true;
		try{
			for(String[] pair:names){
				check(NcrnaFamilyRegistry.displayName(pair[0]).equals(pair[1]), "Display map differs from the selected scientific spelling: "+pair[0]);
				for(int strand=0; strand<2; strand++){for(int partial:new int[]{0, 3}){testEmission(pair[0], pair[1], strand, partial);}}
			}
			final Orf unassigned=fixture(0, 0);check(!unassigned.toGff().contains(",family:"), "An absent optional family attribute must remain absent");
		}finally{GffLine.parseAttributes=oldParse;}
		System.out.println("NCRNA_FAMILY_DISPLAY_NAME_PASS approved_names_both_strands_spaces_following_attributes_partial_internal_ids_unchanged");
	}
	private static void testEmission(String id, String label, int strand, int partial){
		assert(id!=null && label!=null) : "Every fixture binds an independently specified internal key and selected scientific label";
		final Orf orf=fixture(strand, partial);final String withoutFamily=orf.toGff();
		orf.ncrnaFamily=id;final String actual=orf.toGff();
		final int suffix=withoutFamily.indexOf(';');
		final String expected=suffix<0 ? withoutFamily+",family:"+label
			: withoutFamily.substring(0, suffix)+",family:"+label+withoutFamily.substring(suffix);
		check(actual.equals(expected), "Only the family display token may change; preserve coordinates, score, model and partial attributes: "+id);
		check(orf.ncrnaFamily.equals(id) && orf.trnaModel.equals("fixture_model") && orf.ncrnaPartial==partial,
			"Formatting must not mutate internal family, model or endpoint provenance");
		for(String row:new String[]{actual, actual+",following:value"}){
			final GffLine parsed=new GffLine(row.getBytes(StandardCharsets.US_ASCII));
			final String attributes=row.substring(row.lastIndexOf('\t')+1);
			check(attributes.equals(parsed.attributes) && family(parsed.attributes).equals(label),
				"GFF column-nine parsing must preserve the complete scientific name through spaces and later attributes: "+label);
			final GffLine roundTrip=new GffLine(parsed.toString().getBytes(StandardCharsets.US_ASCII));
			check(attributes.equals(roundTrip.attributes) && family(roundTrip.attributes).equals(label)
				&& roundTrip.start==orf.start+1 && roundTrip.stop==orf.stop+1,
				"GffLine serialization and re-parsing must retain coordinates and the entire attributes field, including partial/following fields: "+label);
		}
	}
	private static Orf fixture(int strand, int partial){
		assert((strand==0 || strand==1) && (partial==0 || partial==3)) : "Both strands and complete/physical-edge cases are deliberate fixture dimensions";
		final byte[] bases=new byte[128];Arrays.fill(bases, (byte)'A');
		final Orf orf=new Orf("fixture", 0, 127, strand, 0, bases, false, ProkObject.RNA);
		orf.trnaModel="fixture_model";orf.ncrnaPartial=partial;return orf;
	}
	private static String family(String attributes){
		final String tag=",family:";final int at=attributes.indexOf(tag);
		check(at>=0 && attributes.indexOf(tag, at+tag.length())<0, "Exactly one family attribute must survive serialization");
		final int start=at+tag.length();int stop=start;
		while(stop<attributes.length() && attributes.charAt(stop)!=',' && attributes.charAt(stop)!=';'){stop++;}
		return attributes.substring(start, stop);
	}
	private static void check(boolean condition, String why){if(!condition){throw new AssertionError(why);}}
}
