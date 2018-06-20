#ifndef TRIANGLEGRADIENT
#define TRIANGLEGRADIENT


class TriangleGradient{

    enum TriBitPair {
        ARROWS_EMPTY = 0,
        V0_E01		= 1 << 0
      , V0_E02		= 1 << 1
      , V1_E10		= 1 << 2
      , V1_E12		= 1 << 3
      , V2_E20		= 1 << 4
      , V2_E21		= 1 << 5

      , E01_F012	= 1 << 6
      , E02_F021	= 1 << 7
      , E12_F120	= 1 << 8
    };

private:
    TriBitPair arrows;

public:

    inline TriangleGradient() : arrows(ARROWS_EMPTY) { }

    inline TriangleGradient(unsigned int tetra_gradient)
    {
        arrows = (TriBitPair)(tetra_gradient);
    }

    inline TriBitPair getCode() const{
        return arrows;
    }

    inline bitset<3> getCode(SBitSet s){
        return bitset<3>(s);
    }

    inline bool operator<(TriangleGradient const& c) const{
        return arrows < c.getCode();
    }

    // flag is one of the named entities of the enum
    inline bool testArrow(TriBitPair bit )	const   {  return (arrows & bit) > 0; }                     //  return 0 if flag is not set and 1 if flag is set
    inline void setArrow(TriBitPair bit )           {  arrows = (TriBitPair)(arrows | bit); }  //  sets bit in position corresponding to flag
    inline void clearArrow(TriBitPair bit )         {  arrows = (TriBitPair)(arrows & !bit); }  //  remove bit in position corresponding to flag


    //return the other vertex forming the edge with which (index) is paired. -1 if unpaired
    inline SBitSet get_vertex_pair(SBitSet vertex){
        // use ternary operator to simplify logic here
        switch(vertex){
        case V0:
            return testArrow(V0_E01 )	? E01
                 : testArrow(V0_E02) 	? E02
                 : EMPTY;

        case V1:
            return testArrow(V1_E10 )	? E01
                 : testArrow(V1_E12)	? E12
                 : EMPTY;

        case V2:
            return testArrow(V2_E20 )	? E02
                 : testArrow(V2_E21) 	? E12
                 : EMPTY;

        default:
            cerr << "ERROR, DEFAULT VE" << endl;
            return EMPTY;
        }
    }

    inline bool is_vertex_unpaired(SBitSet vertex){
        return get_vertex_pair(vertex) == EMPTY;
    }


    inline SBitSet get_edge_pair(SBitSet edge){


        switch( edge )
        {
        case E01:
            return testArrow(V0_E01)	? V0
                 : testArrow(V1_E10)    ? V1
                 : testArrow(E01_F012)	? F012
                 : EMPTY;

        case E02:
            return testArrow(V0_E02)	? V0
                 : testArrow(V2_E20)	? V2
                 : testArrow(E02_F021)	? F012
                 : EMPTY;


        case E12:
            return testArrow(V1_E12)			? V1
                 : testArrow(V2_E21)            ? V2
                 : testArrow(E12_F120)			? F012
                 : EMPTY;


        default:
                cerr << "ERROR, DEFAULT EF" << endl;
                return EMPTY;
        }
    }

    inline bool is_edge_unpaired(SBitSet edge){
        return get_edge_pair(edge) == EMPTY;
    }

    inline SBitSet get_triangle_pair(){
        return testArrow(E01_F012)  ? E01
             : testArrow(E02_F021)  ? E02
             : testArrow(E12_F120)  ? E12
             : EMPTY;
    }

    inline bool is_triangle_unpaired(){
        return get_triangle_pair() == EMPTY;
    }

    inline void setVE(SBitSet v, SBitSet e){
        //cout << v << " " << e << endl;
        switch(v){
            case V0:
                switch(e)
                {
                    case E01:	setArrow(V0_E01); break;
                    case E02:	setArrow(V0_E02); break;
                    default: cerr<< "ERROR NEL SET v0"<< endl; break;
                }
                break;
            case V1:
                switch(e)
                {
                    case E01:	setArrow(V1_E10); break;
                    case E12:	setArrow(V1_E12); break;
                    default: cerr<< "ERROR NEL SET v1"<< endl; break;
                }
                break;
            case V2:
                switch(e)
                {
                    case E02:	setArrow(V2_E20); break;
                    case E12:	setArrow(V2_E21); break;
                    default: cerr<< "ERROR NEL SET v2"<< endl; break;
                }
                break;
            default:
                cerr << "ERROR NEL SET VE" << endl;
                break;
        }

    }

    inline void setEF(SBitSet e){
        //cout << e << endl;
        switch(e){
            case E01:
                setArrow(E01_F012); break;
            case E02:
                setArrow(E02_F021); break;
            case E12:
                setArrow(E12_F120); break;
            default:
                cerr << "ERROR NEL SET EF" << endl;
                break;
        }

    }


    ///clear arrows

    inline void freeVE(SBitSet v, SBitSet e){
        //cout << v << " " << e << endl;
        switch(v){
            case V0:
                switch(e)
                {
                    case E01:	clearArrow(V0_E01); break;
                    case E02:	clearArrow(V0_E02); break;
                    default: cerr<< "ERROR NEL SET v0"<< endl; break;
                }
                break;
            case V1:
                switch(e)
                {
                    case E01:	clearArrow(V1_E10); break;
                    case E12:	clearArrow(V1_E12); break;
                    default: cerr<< "ERROR NEL SET v1"<< endl; break;
                }
                break;
            case V2:
                switch(e)
                {
                    case E02:	clearArrow(V2_E20); break;
                    case E12:	clearArrow(V2_E21); break;
                    default: cerr<< "ERROR NEL SET v2"<< endl; break;
                }
                break;
            default:
                cerr << "ERROR NEL SET VE" << endl;
                break;
        }

    }

    inline void freeEF(SBitSet e){
        //cout << e << endl;
        switch(e){
            case E01:
                clearArrow(E01_F012); break;
            case E02:
                clearArrow(E02_F021); break;
            case E12:
                clearArrow(E12_F120); break;
            default:
                cerr << "ERROR NEL SET EF" << endl;
                break;
        }

    }

};

#endif // TRIANGLEGRADIENT

