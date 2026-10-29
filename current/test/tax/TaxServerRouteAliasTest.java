package test.tax;

import java.lang.reflect.Method;
import java.util.HashMap;

import tax.TaxServer;

/** Finite oracle for TaxServer URL route aliases.
 * @author Brian Bushnell, Sucrose */
public final class TaxServerRouteAliasTest{

	/** Runs baseline or fixed expectations for TaxServer route aliases.
	 * @param args Optional expect=baseline or expect=fixed
	 * @throws Exception If reflective map construction fails */
	public static void main(String[] args) throws Exception{
		boolean baseline=false;
		for(String arg : args){
			if(arg.equals("expect=baseline")){baseline=true;}
			else if(arg.equals("expect=fixed")){baseline=false;}
			else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}

		final HashMap<String, Integer> typeMap=typeMap();
		requireAlias(typeMap, "pt_img", EXPECTED_PT_IMG);
		requireAlias(typeMap, "sc_taxid", EXPECTED_SC_TAXID);
		requireAlias(typeMap, "sc_silvaheader", EXPECTED_SC_SILVAHEADER);

		final Integer scImg=typeMap.get("sc_img");
		if(baseline){
			check(scImg==null, "Baseline should not map sc_img; found "+scImg);
		}else{
			check(scImg!=null, "Fixed TaxServer should map sc_img.");
			requireRoute("sc_img", scImg.intValue(), EXPECTED_SC_IMG);
		}

		System.out.println("TAXSERVER_ROUTE_ALIAS_TEST_PASS expect="+(baseline ? "baseline" : "fixed"));
	}

	/** Prevents construction of this standalone fixture. */
	private TaxServerRouteAliasTest(){}

	/** Builds the private TaxServer route-alias map without constructing a server.
	 * @return TaxServer.makeTypeMap() result
	 * @throws Exception If reflection cannot access or run makeTypeMap */
	@SuppressWarnings("unchecked")
	private static HashMap<String, Integer> typeMap() throws Exception{
		final Method method=TaxServer.class.getDeclaredMethod("makeTypeMap");
		method.setAccessible(true);
		return (HashMap<String, Integer>)method.invoke(null);
	}

	/** Requires an alias to exist and encode to the expected route value.
	 * @param typeMap TaxServer alias map
	 * @param alias route alias to inspect
	 * @param expectedRaw expected encoded route value */
	private static void requireAlias(final HashMap<String, Integer> typeMap,
			final String alias, final int expectedRaw){
		final Integer actual=typeMap.get(alias);
		check(actual!=null, "Missing route alias "+alias);
		requireRoute(alias, actual.intValue(), expectedRaw);
	}

	/** Requires a route value to match its independently expected encoded value.
	 * @param alias route alias being checked
	 * @param actual encoded route value from TaxServer
	 * @param expectedRaw expected encoded route value */
	private static void requireRoute(final String alias, final int actual, final int expectedRaw){
		check(actual==expectedRaw,
				alias+" raw route should be "+expectedRaw+"; found "+actual);
	}

	/** Test checks stay active under -da.
	 * @param condition condition required for test success
	 * @param message diagnostic message for failure */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Plaintext IMG route: TaxServer.PT_BIT(16)+TaxServer.IMG(6). */
	private static final int EXPECTED_PT_IMG=22;
	/** Semicolon TaxID route: TaxServer.SC_BIT(32)+TaxServer.TAXID(3). */
	private static final int EXPECTED_SC_TAXID=35;
	/** Semicolon IMG route: TaxServer.SC_BIT(32)+TaxServer.IMG(6). */
	private static final int EXPECTED_SC_IMG=38;
	/** Semicolon SILVA header route: TaxServer.SC_BIT(32)+TaxServer.SILVAHEADER(7). */
	private static final int EXPECTED_SC_SILVAHEADER=39;
}
