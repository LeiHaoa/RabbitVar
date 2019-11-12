   #ifndef _PATTERNS_H
   #define _PATTERNS_H
   
	#define MC_Z_NUM_S_ANY_NUM_S  "\\d+S\\S*\\d+S"

    //Variation patterns
    #define BEGIN_DIGITS  "^(\\d+)"
    #define UP_NUMBER_END "\\^(\\d+)$"
    #define BEGIN_MINUS_NUMBER_ANY "^-\\d+(.*)"
    #define BEGIN_MINUS_NUMBER_CARET "^-\\d+\\^"
    //#define BEGIN_MINUS_NUMBER "^-(\\d+)"
	#define BEGIN_MINUS_NUMBER "^-(\\d+).*"
    #define MINUS_NUM_NUM "-\\d\\d"
    #define HASH_GROUP_CARET_GROUP "#(.+)\\^(.+)"

    //Sclip patterns
    #define B_A7 "^.AAAAAAA.*"
    #define B_T7 "^.TTTTTTT.*"

    //ATGC patterns
    #define CARET_ATGNC "\\^([ATGNC]+).*"
    #define CARET_ATGC_END  "\\^([ATGC]+)$"
    #define AMP_ATGC  ".*&([ATGC]+).*"
    #define BEGIN_PLUS_ATGC "^\\+([ATGC]+).*"
    #define HASH_ATGC "#([ATGC]+).*"
    #define ATGSs_AMP_ATGSs_END "(\\+[ATGC]+)&[ATGC]+$"
    #define MINUS_NUMBER_AMP_ATGCs_END  "(-\\d+)&[ATGC]+$"
    #define MINUS_NUMBER_ATGNC_SV_ATGNC_END "^-\\d+\\^([ATGNC]+)<...\\d+>([ATGNC]+)$"
    #define BEGIN_ATGC_END "^[ATGC]+$"

    //SV patterns
    #define DUP_NUM "<dup(\\d+)"
    #define DUP_NUM_ATGC "<dup(\\d+)>([ATGC]+)$"
    #define INV_NUM "<inv(\\d+)"
    #define SOME_SV_NUMBERS "<(...)\\d+>"
    #define ANY_SV "<(...)>"

    //File and columns patterns
    #define SAMPLE_PATTERN "([^\\/\\._]+).sorted[^\\/]*.bam"
    #define SAMPLE_PATTERN2 "([^\\/]+)[_\\.][^\\/]*bam"
    #define INTEGER_ONLY "^\\d+$"

    //CIGAR patterns
    /**
     * Regexp finds number followed by S followed by number and I or D at the start of string
     */
    #define BEGIN_NUMBER_S_NUMBER_IorD "^(\\d+)S(\\d+)([ID])"
    /**
     * Regexp finds number followed by I or D followed by number followed by S followed by end of string
     */
    #define NUMBER_IorD_NUMBER_S_END "(\\d+)([ID])(\\d+)S$"
    /**
     * Regexp finds number
     * followed by S followed by number followed by M followed by number followed by I or D at the start of string
     */
    #define BEGIN_NUMBER_S_NUMBER_M_NUMBER_IorD "^(\\d+)S(\\d+)M(\\d+)([ID])"
    /**
     * Regexp finds number
     * followed by I or D followed by number followed by M followed by number followed by S followed by end of string
     */
    #define NUMBER_IorD_NUMBER_M_NUMBER_S_END "(\\d+)([ID])(\\d+)M(\\d+)S$"
    /**
     * Regexp finds digit-M-number-(I or D)-number-M
     */
    #define BEGIN_DIGIT_M_NUMBER_IorD_NUMBER_M "^(\\d)M(\\d+)([ID])(\\d+)M"

    #define BEGIN_DIGIT_M_NUMBER_IorD_NUMBER_M_ "^\\dM\\d+[ID]\\d+M"
    /**
     * Regexp replaces number-(I or D)-number-M
     */
    #define NUMBER_IorD_DIGIT_M_END "(\\d+)([ID])(\\d)M$"
    #define NUMBER_IorD_NUMBER_M_END "(\\d+)([ID])(\\d+)M$"
    /**
     * Regexp finds number-M-number-D-digit-M-number-I-digit-M-number-D. ^(.*?) captures everything before the
     * M-D-M-I-M-D complex
     */
    #define D_M_D_DD_M_D_I_D_M_D_DD "^(.*?)(\\d+)M(\\d+)D(\\d+)M(\\d+)I(\\d+)M(\\d+)D(\\d+)M"
    #define D_M_D_DD_M_D_I_D_M_D_DD_prim "(\\d+)M(\\d+)D(\\d+)M(\\d+)I(\\d+)M(\\d+)D(\\d+)M"
    //------------[haoz:] ??? -------------
    #define threeDeletionsPattern "^(.*?)(\\d+)M(\\d+)D(\\d+)M(\\d+)D(\\d+)M(\\d+)D(\\d+)" + "M"
    #define threeIndelsPattern "^(.*?)(\\d+)M(\\d+)([DI])(\\d+)M(\\d+)([DI])(\\d+)M(\\d+)([DI])(\\d+)M"
    #define DIGM_D_DI_DIGM_D_DI_DIGM_DI_DIGM "\\d+M\\d+[DI]\\d+M\\d+[DI]\\d+M\\d+[DI]\\d+M"
    #define DM_DD_DM_DD_DM_DD_DM "\\d+M\\d+D\\d+M\\d+D\\d+M\\d+D\\d+M"
    /**
     * Regexp finds number-D-digit-M-number-D and optionally number-I
     */
    #define DIG_D_DIG_M_DIG_DI_DIGI "(\\d+)D(\\d+)M(\\d+)([DI])(\\d+I)?"
    /**
     * Regexp finds number-I-digit-M-number and optionally number-I
     */
    #define DIG_I_DIG_M_DIG_DI_DIGI "(\\d+)I(\\d+)M(\\d+)([DI])(\\d+I)?"
    /**
     * Regexp finds not_digit-number-I-number-M-number-D and optionally number-I
     */
    #define NOTDIG_DIG_I_DIG_M_DIG_DI_DIGI "(\\D)(\\d+)I(\\d+)M(\\d+)([DI])(\\d+I)?"
    /**
     * Regexp finds number-D-number-D
     */
    #define DIG_D_DIG_D "(\\d+)D(\\d+)D"
    /**
     * Regexp finds number-I-number-I
     */
    #define DIG_I_DIG_I "(\\d+)I(\\d+)I"
    #define BEGIN_ANY_DIG_M_END "^(.*?)(\\d+)M$"

    #define DIG_M_END "\\d+M$"
    #define BEGIN_DIG_M "^(\\d+)M"
    /**
     * Regexp finds number-S-number-M at start of CIGAR string
     */
    #define DIG_S_DIG_M "^(\\d+)S(\\d+)M"
    #define DIG_M_DIG_S_END "\\d+M\\d+S$"
    /**
     * Regexp finds number-M-number-S at end of CIGAR string ^(.*?) captures everything before M-S complex
     */
    #define ANY_NUMBER_M_NUMBER_S_END "^(.*?)(\\d+)M(\\d+)S$"
    /**
     * Regexp finds number followed by D at the start of string
     */
    #define BEGIN_NUMBER_D "^(\\d+)D"
    #define END_NUMBER_D "(\\d+)D$"
    /**
     * Regexp finds number followed by I at the start of string
     */
    #define BEGIN_NUMBER_I "^(\\d+)I"
    /**
     * Regexp finds number followed by I at the end of string
     */
    #define END_NUMBER_I "(\\d+)I$"

    /**
     * Regexp finds numbers followed by M (matched) or D (deleted) in CIGAR string
     */
    #define ALIGNED_LENGTH_MND "(\\d+)[MND]"

    /**
     * The total aligned length, excluding soft-clipped bases and insertions
     */
    #define ALIGNED_LENGTH_MD "(\\d+)[MD=X]"

    #define SOFT_CLIPPED "(\\d+)[MIS]"
    #define SA_CIGAR_D_S_5clip "^\\d\\d+S"
    #define SA_CIGAR_D_S_3clip "\\d\\dS$"
    #define BEGIN_dig_dig_S_ANY_dig_dig_S_END "^\\d\\dS.*\\d\\dS$"
    #define BEGIN_NUM_S_OR_BEGIN_NUM_H "^(\\d+)S|^\\d+H"
    #define END_NUM_S_OR_NUM_H "(\\d+)S$|H$"

#endif
