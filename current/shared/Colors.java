package shared;

/**
 * Utility class for ANSI terminal color codes and text formatting.
 * Provides methods for applying colors, underlines, and bold formatting to text output.
 * Automatically detects Windows compatibility and adjusts output accordingly.
 * @author Brian Bushnell
 */
public class Colors {
	
	/**
	 * Test method demonstrating color formatting capabilities.
	 * Prints sample colored text to demonstrate red, blue, green, and yellow colors.
	 * @param args Command-line arguments (not used)
	 */
	public static void main(String[] args){
		System.out.println(format("Red", RED, false));
		System.out.println(format("Blue", BLUE, false));
		System.out.println(format("Green", GREEN, true));
		System.out.println(format("Yellow", YELLOW, true));
	}
	
	/**
	 * Formats a string with ANSI color codes and optional underline.
	 * Applies the specified color and underline formatting, then resets to default.
	 *
	 * @param s The string to format
	 * @param color ANSI color code to apply
	 * @param underline true to add underline formatting, false otherwise
	 * @return Formatted string with color codes, or original string if colors disabled
	 */
	public static String format(String s, String color, boolean underline){
//		return color+(bold ? "[1m":"")+(bold ? "[4m":"")+s+RESET;
		return color+(underline ? UNDERLINE : "")+s+RESET;
	}
	
	/**
	 * Formats a byte array as colored text with optional underline.
	 * Converts byte array to string, applies formatting, then converts back to bytes.
	 *
	 * @param s Byte array containing text to format
	 * @param color ANSI color code to apply
	 * @param underline true to add underline formatting, false otherwise
	 * @return Formatted text as byte array
	 */
	public static byte[] format(byte[] s, String color, boolean underline){
		return format(new String(s), color, underline).getBytes();
	}
	
	/**
	 * Creates an array containing all available colors (bright and standard).
	 * Includes bright colors first, then standard colors for maximum variety.
	 * @return Array of ANSI color codes in priority order
	 */
	public static String[] makeColorArray(){
		return new String[] {
				BRIGHT_RED, BRIGHT_GREEN, BRIGHT_YELLOW, BRIGHT_BLUE, BRIGHT_PURPLE, BRIGHT_CYAN, RED, GREEN, YELLOW, PURPLE, CYAN
		};
	}
	
	/**
	 * Creates an array containing only standard (dark) color codes.
	 * Excludes bright variants for subdued color output.
	 * @return Array of standard ANSI color codes
	 */
	public static String[] makeDarkArray(){
		return new String[] {
				RED, GREEN, YELLOW, BLUE, PURPLE, CYAN,
		};
	}
	
	/**
	 * Creates an array containing only bright color codes.
	 * Provides high-contrast colors for emphasized output.
	 * @return Array of bright ANSI color codes
	 */
	public static String[] makeBrightArray(){
		return new String[] {
				BRIGHT_RED, BRIGHT_GREEN, BRIGHT_YELLOW, BRIGHT_BLUE, BRIGHT_PURPLE, BRIGHT_CYAN
		};
	}

	/** Flag indicating whether the system is running Windows */
	public static boolean windows = Shared.WINDOWS;
	/** Flag to skip color output (true on Windows due to limited ANSI support) */
	public static boolean skip = windows;
	/**
	 * ANSI escape sequence prefix (literal text on Windows, actual escape on Unix)
	 */
	public static String esc = windows ? "<ESC>" : "\u001B"; //Windows only supports colors with Win10+

	/** ANSI code to reset all text formatting to default */
	public static String RESET = skip ? "" : esc+"[0m";
	/** ANSI code to enable underline text formatting */
	public static String UNDERLINE = skip ? "" : esc+"[4m";
	/** ANSI code to enable bold text formatting */
	public static String BOLD = skip ? "" : esc+"[1m";
	
	/** ANSI color code for black text */
	public static String BLACK = skip ? "" : esc+"[30m";
	/** ANSI color code for red text */
	public static String RED = skip ? "" : esc+"[31m";
	/** ANSI color code for green text */
	public static String GREEN = skip ? "" : esc+"[32m";
	/** ANSI color code for yellow text */
	public static String YELLOW = skip ? "" : esc+"[33m";
	/** ANSI color code for blue text */
	public static String BLUE = skip ? "" : esc+"[34m";
	/** ANSI color code for purple text */
	public static String PURPLE = skip ? "" : esc+"[35m";
	/** ANSI color code for cyan text */
	public static String CYAN = skip ? "" : esc+"[36m";
	/** ANSI color code for white text */
	public static String WHITE = skip ? "" : esc+"[37m";
	
	/** ANSI color code for bright black (dark gray) text */
	public static String BRIGHT_BLACK = skip ? "" : esc+"[30;1m";
	/** ANSI color code for bright red text */
	public static String BRIGHT_RED = skip ? "" : esc+"[31;1m";
	/** ANSI color code for bright green text */
	public static String BRIGHT_GREEN = skip ? "" : esc+"[32;1m";
	/** ANSI color code for bright yellow text */
	public static String BRIGHT_YELLOW = skip ? "" : esc+"[33;1m";
	/** ANSI color code for bright blue text */
	public static String BRIGHT_BLUE = skip ? "" : esc+"[34;1m";
	/** ANSI color code for bright purple text */
	public static String BRIGHT_PURPLE = skip ? "" : esc+"[35;1m";
	/** ANSI color code for bright cyan text */
	public static String BRIGHT_CYAN = skip ? "" : esc+"[36;1m";
	/** ANSI color code for bright white text */
	public static String BRIGHT_WHITE = skip ? "" : esc+"[37;1m";

	/** Pre-created array containing all available colors (bright and standard) */
	public static String[] colorArray=makeColorArray();
	/** Pre-created array containing only standard (dark) color codes */
	public static String[] darkArray=makeDarkArray();
	/** Pre-created array containing only bright color codes */
	public static String[] BrightArray=makeBrightArray();
	
}
