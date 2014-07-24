import os, sys, re
import argparse
from xml.etree.ElementTree import iterparse

PROTEIN     = 'protein'
INTERACTION = 'interaction'
SIGNATURE   = 'signature'
POSITIVE    = 'positive'
NEGATIVE    = 'negative'
LOOP        = 'loop'
DOMAIN      = 'domain'
MAPPING     = 'mapping'
ALIGNMENT   = 'alignment'
RFRESULT    = 'Random Forest result'
ALL         = 'all'

class ILXMLParser(object): 
    '''
    STRUCTURE OF XML FILE 
    
    ---------------------------------------------------------------------------------------------------------
    Level                                    XMLtags                                 Object/Sub-Object/Buffer
    ---------------------------------------------------------------------------------------------------------

    R: --XML.................................<xml/>..................................ILXMLParser (Buffer)
    0:   +--Protein..........................<protein/>...................................+ILXMLProtein (Object)
         |  +--Name..........................<name/>                                      |     +--name
         |  +--MD5...........................<MD5/>                                       |     +--MD5
         |  +--Sequence......................<sequence/>                                  |     +--Sequence
         |  +--{Alignments}..................<blastSearch><search_type/></blastSearch>    |     |
         |  |   +--AlignmentType.............<ID/>                                        |     |
         |  |   +--AlignmentIteration........<search_iteration/>                          |     |
         |  |   +--..........................<alignments/>                                |     +--{Alignments}
    1:   |  |   +--Alignment.................<align/>.....................................|.....|...+--ILXMLAlignment (Sub-Object)
         |  |      |                                                                      |     |           +--ali_type (<ID>/>)
         |  |      |                                                                      |     |           +--iter_no (<search_iteration/>)
         |  |      +TargetID.................<target_id/>                                 |     |           +--target_id
         |  |      +eValue...................<e-value/>                                   |     |           +--eValue
         |  |      +queryStart...............<query_start/>                               |     |           +--q_start
         |  |      +queryEnd.................<query_end/>                                 |     |           +--q_end
         |  |      +querySequence............<query_seq/>                                 |     |           +--q_seq
         |  |      +hitStart.................<hit_start/>                                 |     |           +--h_start
         |  |      +hitEnd...................<hit_end/>                                   |     |           +--h_end
         |  |      +hitSequence..............<hit_seq/>                                   |     |           +--h_seq
         |  |                                                                             |     |
         |  +--{DomainMappings}..............<domainMapping><alignments/></domainMapping> |     +--{DomainMappings}
         |  |   +--DomainMapping.............<align/>.....................................|.....|...+--ILXMLDomainMapping (Sub-Object)
         |  |      +--targetID...............<target_id/>                                 |     |           +--target_id
         |  |      +--eValue.................<e-value/>                                   |     |           +--eValue
         |  |      +--queryStart.............<query_start/>                               |     |           +--q_start
         |  |      +--queryEnd...............<query_end/>                                 |     |           +--q_end
         |  |      +--hitStart...............<hit_start/>                                 |     |           +--h_start
         |  |      +--hitEnd.................<hit_end/>                                   |     |           +--h_end
         |  |                                                                             |     |
         |  +--{Features}....................<features/>                                  |     +--{Features}
    1:   |      +--Feature...................<feature/>...................................|.........+--ILXMLFeature ([ILXMLLoop, ILXMLDomain]) (Sub-Object)
         |         +--type...................<type/>                                      |                 +--Type ([Loop, Domain])
         |         +--name...................<loop_subClass/>|<SCOP_code/>                |                 +--([subClass, SCOP_code])
         |         +--structureMatch.........<loopMatch/>    |<SCOP_match/>               |                 +--([loopMatch, SCOP_match])
         |         +--alignmentIni...........<align_ini/>                                 |                 +--ali_ini
         |         +--alignmentEnd...........<align_end/>                                 |                 +--ali_end
         |         +--queryFeatureSeq........<query_feature_ali_seq/>                     |                 +--QAliSeq
         |         +--targetFeatureSeq.......<target_feature_ali_seq/>                    |                 +--HAliSeq
         |         +--{domainAssignations}...<SCOP_assignations/>|<domain_mapping/>       |                 +--{mappings} (only ILXMLLoop)
    2:   |             +--domainAssignation..<SCOP_assignation/> |<domain_mapping/>.......|.....................+ILXMLDomainAssignation (Buffer)
         |                +--type............<SCOP_assignation/> |<domain_type/>          |                           +--Type
         |                +--ID..............<SCOP_id/>          |<domain_id/>            |                           +--ID 
         |                +--ini.............<SCOP_ini/>         |<domain_ini/>           |                           +--ini
         |                +--end.............<SCOP_end/>         |<domain_end/>           |                           +--end
         |                                                                                |
    0:   +--Interaction......................<interaction/>...............................+ILXMLInteraction (Object)
            +--Protein1ID....................<P1ID/>                                            +--interactor1_ID (index in ILXMLParser.proteins)
            +--Protein2ID....................<P2ID/>                                            +--interactor2_ID (index in ILXMLParser.proteins)
            +--Splus.........................<Splus/>                                           +--Splus
            +--Sminus........................<Sminus/>                                          +--Sminus
            +--LSR...........................<LSR/>                                             +--LSR
            +--LpVR..........................<LpVR/>                                            +--LpVR
            +--{Random Forest Results}.......<RF_results/>                                      +--{RF_results}
    1:      |   +--Random Forest Result......<RF_result/>.......................................|...+ILXMLInteractionRFResult (Sub-Object)
            |      +--Relative Cost..........<RF_cost/>                                         |         +--cost
            |      +--Prediction.............<RF_prediction/>                                   |         +--prediction <bool>
            |      +--Random Forest Score....<RF_score/>                                        |         +--RF_score
            |      +--{Precisions}...........<inferred_precisions/>                             |         +--{precisions}
    2:      |          +--Precision..........<inferred_precision/>..............................|.............+--ILXMLPrecision (Buffer)
            |             +--Unbalance Ratio.<unbalance_ratio/>                                 |                     +--UR
            |             +--Precision value.<precision/>                                       |                     +--precision
            |             +--Precision error.<error/>                                           |                     +--error
            +--{Interaction Signatures}......<interaction_signatures/>                          +--{signatures}
    1:          +--Interaction Signature.....<interaction_signature/>...............................+ILXMLInteractionSignature (Sub-Object)
                   +--Signature Type.........<IS_type/>                                                   +refMatrix
                   |                                                                                      +elementsType ([LOOP, DOMAIN])
                   +--Protein1 Sign. ID......<P1IS/>                                                      +P1SignID (<tuple>, ref. to self.loops or 
                   |                                                                                      |                           self.domains index)
                   +--Protein2 Sign. ID......<P2IS/>                                                      +P2SignID (<tuple>, ref. to self.loops or 
                   |                                                                                      |                           self.domains index)
                   +--pValue.................<pValue/>                                                    +pValue
    '''


    def __init__(self): 
        self._loops = set()
        self._domains = set()
        self.loops = []
        self.domains = []
        self.proteins = []
        self.interactions = []
        self.active_object = self.reset_active_object()
        self.active_subObj = self.reset_active_subObj()
        self.kwds = self.reset_kwds()
        self.last_output_level = -1
        self.skip_levels = -1

    # SETTERS: CLEAN UP METHODS

    def reset_active_object(self): self.active_object = None
    def reset_active_subObj(self): self.active_subObj = None
    def reset_kwds(self)         : self.kwds          = {}


    # SETTERS ACTIVE OBJECT/SUB-OBJECT CONTROL METHODS

    def set_active_object(self, objType): 
        if   objType == PROTEIN: 
            self.active_object = ILXMLProtein(ILXMLParser=self, ID=len(self.get_proteins()))

        elif objType == INTERACTION: 
            self.active_object = ILXMLInteraction(ILXMLParser=self, ID=len(self.get_interactions()))
            
    def set_active_subObj(self, subObj_generator, parentID, parentType, **kwds): 
        self.active_subObj = subObj_generator(self, parentID, parentType, **kwds)


    # SETTERS: METHODS FOR STORING (REFERENCES TO) OBJECTS/SUB-OBJECTS IN THE PARSER

    def add_protein(self, proteinObj, store=True):
        '''
        o proteinObj: ILXMLProtein object
        o store     : <bool>
        '''
        if store is True: self.proteins.append(proteinObj)
        else            : self.proteins.append(proteinObj.get_name())

    def add_interaction(self, interactionObj, store=True):
        '''
        o proteinObj: ILXMLInteraction object
        o store     : <bool>
        '''
        if store is True: self.interactions.append(interactionObj)
        else            : self.interactions.append( (interactionObj.get_i1name(), interactionObj.get_i2name()) )

    
    def add_feature(self, code, feature_type):
        if   feature_type == LOOP: 
            ref_set  = self._loops
            ref_list = self.get_loops()
        elif feature_type == DOMAIN: 
            ref_set  = self._domains
            ref_list = self.get_domains()

        if not code in ref_set: 
            ref_set.add(code)
            ref_list.append(code)
            
        return ref_list.index(code)


    # PRIVATE SETTERS: METHODS FOR CONTROLLING THE OUTPUT BEHAVIOUR
    def _set_kwds(self, kwds)             : self.kwds              = kwds
    def set_last_output_level(self, level): self.last_output_level = level
    def set_skip_levels(self, level_no)   : self.skip_levels       = level_no


    # GETTERS

    def get_proteins(self)         : return self.proteins
    def get_interactions(self)     : return self.interactions
    def get_loops(self)            : return self.loops
    def get_domains(self)          : return self.domains
    def get_active_object(self)    : return self.active_object
    def get_active_subObj(self)    : return self.active_subObj
    def get_last_output_level(self): return self.last_output_level
    def get_skip_levels(self)      : return self.skip_levels


    # METHODS FOR INITIALIZING OBJECTS/SUB-OBJECTS/BUFFERS
    def _start_protein(self)                   : self.set_active_object( PROTEIN )
    def _start_interaction(self)               : self.set_active_object( INTERACTION )
    def _start_alignment(self, parentProteinID, ali_type, ali_iter): self.set_active_subObj( subObj_generator = ILXMLAlignment, 
                                                                                             parentID         = parentProteinID, 
                                                                                             parentType       = PROTEIN, 
                                                                                             ali_type         = ali_type, 
                                                                                             iter_no          = ali_iter)
    def _start_domain_mapping(self, parentProteinID): self.set_active_subObj( subObj_generator = ILXMLDomainMapping, 
                                                                              parentID         = parentProteinID, 
                                                                              parentType       = PROTEIN)
    def _start_loop_feature(self, parentProteinID): self.set_active_subObj( subObj_generator = ILXMLLoop, 
                                                                            parentID         = parentProteinID, 
                                                                            parentType       = PROTEIN)
    def _start_domain_feature(self, parentProteinID): self.set_active_subObj( subObj_generator = ILXMLDomain, 
                                                                              parentID         = parentProteinID, 
                                                                              parentType       = PROTEIN)

    def _start_RF_result(self, parentInteractionID): self.set_active_subObj( subObj_generator = ILXMLInteractionRFResult, 
                                                                             parentID         = parentInteractionID, 
                                                                             parentType       = INTERACTION)

    def _start_interaction_signature(self, parentInteractionID): self.set_active_subObj( subObj_generator = ILXMLInteractionSignature, 
                                                                                         parentID         = parentInteractionID, 
                                                                                         parentType       = INTERACTION)

    # METHODS FOR ID/HUMAN READABLE CODE CONVERSION
    def get_proteinID(self, protein_code): 
        try: return self.get_proteins().index(protein_code)
        except ValueError:
            protein_codes = [ x.get_name() for x in self.get_proteins() ]
            return protein_codes.index(protein_code)
                              

    def get_featureCode(self, featureID, feature_type): 
        if   feature_type == LOOP  : ref_list = self.get_loops()
        elif feature_type == DOMAIN: ref_list = self.get_domains()
        return ref_list[featureID]

    def get_featureID(self, feature_code, feature_type=None): 
        if feature_type is None: 
            if   feature_code in self._loops()  : feature_type = LOOP
            elif feature_code in self._domains(): feature_type = DOMAIN

        if   feature_type == LOOP  : ref_list = self.get_loops()
        elif feature_type == DOMAIN: ref_list = self.get_domains()
        return ref_list.index(feature_code)
    
    def signatureCode_2_singatureID(self,feature_codes_list):
        if   eval(feature_codes_list)[0] in self._loops  : feature_type = LOOP
        elif eval(feature_codes_list)[0] in self._domains: feature_type = DOMAIN
        else: 
            print(self._domains)
            print(self._loops)
            raise ValueError(eval(feature_codes_list)[0])
        
        return [ feature_type, [ self.get_featureID(feature_code, feature_type) for feature_code in eval(feature_codes_list) ]]


    def signatureID_2_signatureCode(self, feature_type, feature_IDs_list): 
        if   feature_type == LOOP  : ref_list = self.get_loops()
        elif feature_type == DOMAIN: ref_list = self.get_domains()

        return [ ref_list[x] for x in feature_IDs_list ]
    

    # OUTPUT CONTROL METHODS
    def _check_output_error(self, object_instance_checkClass, object_instance_str, output_level, output_object=None): 
        '''
        o object_instance_checkClass: class required for the active object or sub-object (depending on output level)
        o object_instance_str       : <str> string describing the class required for the active object or sub-object (depending on output level)
        o output_level              : <int> [0,1]. 0 for output of objects. 1 for output of sub-objects
        o output_object             : object to output
                                      If None (default) the method checks the:
                                      a) self.get_active_object() if output_level is 0
                                      b) self.get_active_subObj() if output_level is 1 
        '''
        if   output_level == 0: 
            obj_str   = 'object'
            check_obj = self.get_active_object()
        elif output_level == 1: 
            obj_str   = 'sub-object'
            check_obj = self.get_active_subObj()

        if   not output_object is None and not isinstance(output_object, object_instance_checkClass): 
            err_msg = 'provided %s must be an instance of %s class(es). Obtained %s is %s\n' %(obj_str, object_instance_str, obj_str, type(check_obj))
            raise ValueError(err_msg)
        elif output_object is None and not isinstance(check_obj, object_instance_checkClass): 
            err_msg = 'active %s must be an instance of %s class(es). Obtained %s is %s\n' %(obj_str, object_instance_str, obj_str, type(check_obj))
            raise ValueError(err_msg)
        elif not output_object is None: report_object = output_object
        else                          : report_object = check_obj

        return report_object

    
    def _control_skip_output(self, output_level):
        if self.get_skip_levels() == output_level: 
            self.set_skip_levels(output_level-1)
            self.set_last_output_level(output_level) # upon return, this level will be skipped, but still need to be annotated as outputed!
            return True

        return False


    def _output_protein_info(self, report_level, protein_object = None, 
                             output_alignments=True, output_domain_mappings=True, output_protein_features=True, output_domain_assignations=True): 
        if self._control_skip_output(output_level = 0) is True: return

        report_object = self._check_output_error(object_instance_checkClass = ILXMLProtein, 
                                                 object_instance_str        = 'ILXMLProtein', 
                                                 output_level               = 0, 
                                                 output_object              = protein_object)

        yield( self.custom_protein_output( report_object, **self.kwds ) )
        if report_level < 1: 
            if output_alignments is True: 
                for alignment in report_object.get_alignments()   : yield( self.custom_alignment_output(alignment, **self.kwds) )

            if output_domain_mappings is True: 
                for domain_mapping in report_object.get_mappings(): yield( self.custom_domain_mapping_output(domain_mapping, **self.kwds) )

            if output_protein_features is True: 
                for domain in report_object.get_domains()         : yield( self.custom_feature_output(domain, DOMAIN, **self.kwds) )
                for loop in report_object.get_loops()             : 
                    yield( self.custom_feature_output(loop, LOOP, **self.kwds) )
                    for domain_assignation in loop.get_mappings(): yield( self.custom_domain_assignation_output(domain_assignation, **self.kwds) )

        self.set_last_output_level(0)
                     

    def _output_protein_alignment_info(self, report_level, alignment_object=None): 
        if self._control_skip_output(output_level = 1) is True: return

        report_object = self._check_output_error(object_instance_checkClass = ILXMLAlignment, 
                                                 object_instance_str        = 'ILXMLAlignment', 
                                                 output_level               = 1, 
                                                 output_object              = alignment_object)

        if self.get_last_output_level() < 1:
            yield( self.custom_protein_output(self.get_active_object(), **self.kwds) )
            self.set_skip_levels(0)
        
        yield( self.custom_alignment_output(report_object, **self.kwds) )
        self.set_last_output_level(1)


    def _output_protein_domain_mapping_info(self, report_level, domain_mapping_object=None, retrieve_results=False):
        if self._control_skip_output(output_level = 1) is True: return

        report_object = self._check_output_error(object_instance_checkClass = ILXMLDomainMapping, 
                                                 object_instance_str        = 'ILXMLDomainMapping', 
                                                 output_level               = 1, 
                                                 output_object              = domain_mapping_object)

        if self.get_last_output_level() < 1: 
            yield( self.custom_protein_output(self.get_active_object(), **self.kwds) )
            self.set_skip_levels(0)

        yield( self.custom_domain_mapping_output(report_object, **self.kwds) )
        self.set_last_output_level(1)


    def _output_protein_features_info(self, report_level, feature_object=None, retrieve_results=False, output_domain_assignations=True):
        if self._control_skip_output(output_level = 1) is True: return

        report_object = self._check_output_error(object_instance_checkClass = ILXMLFeature, 
                                                 object_instance_str        = 'ILXMLLoop or ILXMLDomain', 
                                                 output_level               = 1, 
                                                 output_object              = feature_object)

        if self.get_last_output_level() < 1: 
            yield( self.custom_protein_output(self.get_active_object(), **self.kwds) )
            self.set_skip_levels(0)

        if   isinstance(report_object, ILXMLDomain): yield( self.custom_domain_output(report_object, **self.kwds) )
        elif isinstance(report_object, ILXMLLoop): 
            yield( self.custom_loop_output(report_object, **self.kwds) )
            if report_level < 2: 
                for domain_assignation in report_object.get_mappings(): 
                    yield( self.custom_domain_assignation_output(domain_assignation, **self.kwds) )

        self.set_last_output_level(1)


    def _output_feature_domain_assignation_info(self, report_level, domain_assignation_object, retrieve_results=False): 
        if self._control_skip_output(output_level = 2) is True: return

        if self.get_last_output_level() < 1: 
            yield( self.custom_protein_output(self.get_active_object(), **self.kwds) )
            self.set_skip_levels(0)

        if self.get_last_output_level() < 2: 
            if   isinstance(self.get_active_subObj(), ILXMLLoop)  : yield( self.custom_loop_output(self.get_active_subObj(), **self.kwds) )
            elif isinstance(self.get_active_subOjb(), ILXMLDomain): yield( self.custom_domain_ouptut(self.get_active_subObj(), **self.kwds) )
            self.set_skip_levels(1)
        
        yield( self.custom_RF_precision_output(domain_assignation_object, **self.kwds) )
        self.set_last_output_level(2)


    def _output_interaction_info(self, report_level, interaction_object=None, retrieve_results=False, 
                                 output_interaction_signatures=True, output_RF_results=True, output_RF_precisions=True):     
        if self._control_skip_output(output_level = 0) is True: return

        report_object = self._check_output_error(object_instance_checkClass = ILXMLInteraction, 
                                                 object_instance_str        = 'ILXMLInteraction', 
                                                 output_level               = 0, 
                                                 output_object              = interaction_object)

        yield( self.custom_interaction_output( report_object, **self.kwds ) )

        if report_level < 1:
            if output_interaction_signatures is True: 
                for posSign in report_object.get_positive_signatures(): 
                    yield( self.custom_interaction_signature_output(posSign, POSITIVE, **self.kwds) )
                for negSign in report_object.get_negative_signatures(): 
                    yield( self.custom_interaction_signature_output(negSign, NEGATIVE, **self.kwds) )

            if output_RF_results is True: 
                for RFResult in report_object.get_RFResults(): 
                    yield( self.custom_RF_score_output(RFResult, **self.kwds) )
                    if output_RF_precisions is True: 
                        for precision in RFResult.get_precisions(): yield( self.custom_RF_precision_output(precision, **self.kwds) )

        self.set_last_output_level(0)


    def _output_RF_result_info(self, report_level, RF_result_object=None, retrieve_results=False, output_RF_precisions=True):  
        if self._control_skip_output(output_level = 1) is True: return
    
        report_object = self._check_output_error(object_instance_checkClass = ILXMLInteractionRFResult, 
                                                 object_instance_str        = 'ILXMLInteractionRFResult', 
                                                 output_level               = 1, 
                                                 output_object              = RF_result_object)
        if self.get_last_output_level() < 1: 
            yield( self.custom_interaction_output(self.get_active_object(), **self.kwds) )
            self.set_skip_levels(0)

        yield( self.custom_RF_score_output(report_object, **self.kwds) )

        if report_level < 2: 
            for precision in report_object.get_precisions(): yield( self.custom_RF_precision_output(precision, **self.kwds) )

        self.set_last_output_level(1)


    def _output_RF_precision_info(self, report_level, RF_precision_object, retrieve_results=False): 
        if self._control_skip_output(output_level = 2) is True: return
        
        if self.get_last_output_level() < 1: 
            yield( self.custom_interaction_output(self.get_active_object(), **self.kwds) )
            self.set_skip_levels(0)

        if self.get_last_output_level() < 2: 
            yield( self.custom_RF_score_output(self.get_active_subObj(), **self.kwds) )
            self.set_skip_levels(1)

        yield( self.custom_RF_precision_output(RF_precision_object, **self.kwds) )
        self.set_last_output_level(2)


    def _output_interaction_signature_info(self, report_level, interaction_signature_object=None, retrieve_results=False): 
        if self._control_skip_output(output_level = 1) is True: return
 
        report_object = self._check_output_error(object_instance_checkClass = ILXMLInteractionSignature, 
                                                 object_instance_str        = 'ILXMLInteractionSignature', 
                                                 output_level               = 1, 
                                                 output_object              = interaction_signature_object)

        if self.get_last_output_level() < 1: 
            yield( self.custom_interaction_output(self.get_active_object(), **self.kwds) )
            self.set_skip_levels(0)
        
        yield( self.custom_interaction_signature_output(report_object, report_object.get_refMatrix(), **self.kwds) )
        self.set_last_output_level(1)


    # CUSTOM OUTPUT METHODS (should be overwritten by user for output style modification) 
    def custom_protein_output(self, protein_object, **kwds)                                    : return repr(protein_object)    
    def custom_alignment_output(self, alignment_subObj, **kwds)                                : return repr(alignment_subObj)
    def custom_domain_mapping_output(self, domain_mapping_subObj, **kwds)                      : return repr(domain_mapping_subObj)
    def custom_feature_output(self, feature_subObj, feature_type, **kwds)                      : return repr(feature_subObj)
    def custom_domain_assignation_output(self, domain_assign_buffer, **kwds)                   : return repr(domain_assign_buffer)
    def custom_interaction_output(self, interaction_object, **kwds)                            : return repr(interaction_object)
    def custom_interaction_signature_output(self, int_signature_subObj, signature_type, **kwds): return repr(int_signature_subObj)
    def custom_RF_score_output(self, RF_score_subObj, **kwds)                                  : return repr(RF_score_subObj)
    def custom_RF_precision_output(self, RF_precision_buffer, **kwds)                          : return repr(RF_precision_buffer)

    def custom_domain_output(self, domain_subObj, **kwds): return self.custom_feature_output(domain_subObj, DOMAIN, **kwds)
    def custom_loop_output(self, loop_subObj, **kwds)    : return self.custom_feature_output(loop_subObj, LOOP, **kwds)


    # PARSER METHODS
    def results_parser(self, xml_file, report_level=0, skip_alignment_type=[], output_proteins=True, 
                       output_alignments=True, output_domain_mappings=True, output_protein_features=True, output_domain_assignations=True, 
                       output_interactions=True, output_interaction_signatures=True, output_RF_results=True, output_RF_precisions=True, **kwds):
        
        '''
        o xml_file                     : <str>  iLoops XML results file to be parsed
        o report_level                 : <int>  contained in set([-1,0,1,2])
                                         If -1: Only report after parsing the full XML file. Stores all information of the XML document into memmory. 
                                         If  0: Only report after reading a complete protein or interaction. 
                                         If  1: Report after reading a level 1 object 
                                                (Alignment, DomainMapping, Feature for proteins)
                                                (InteractionSignature, RFResult for interactions)
                                         If  2: Report after reading a level 2 object
        o skip_alignment_type          : <list> containing all the types of alignments to be skiped during parse. 
                                         Use [ALL,] to skip all (do not parse alignments information)
        o output_proteins              : <bool> If True, info about proteins is outputed through self.custom_protein_output() method
        o output_alignments            : <bool> If True, info about alignments is outputed through self.custom_alignment_output() method
        o output_domain_mappings       : <bool> If True, info about domain mappings is outputed through self.custom_domain_mapping_output() method
                                                Domain mappins are domains mapped to the whole protein. 
        o output_protein_features      : <bool> If True, info about features (loops, domains) is outputed through:
                                                         Domains: self.custom_loop_output()
                                                         Loops  : self.custom_domain_output()
                                                A generic self.custom_fature_output() is also available. 
        o output_domain_assignations   : <bool> If True, info about domain assignations is outputed through self.custom_domain_assignation_output() method. 
                                                Domain assignations are domains assigned (mapped) to loop features
        o output_interactions          : <bool> If True, info about interactions is outputed through self.custom_interaction_output() method
        o output_interaction_signatures: <bool> If True, info about int. signatures is outputed through self.custom_interaction_signature_output() method.
        o output_RF_results            : <bool> If True, info about random forest results is outputed trhough self.custom_RF_score_output() method
        o output_RF_precisions         : <bool> If True, info about unabalance ratios and inferred precisions of the random forest results is outputed
                                                         through self.custom_RF_precision_output method()
        '''

        # SET STORE CONDITIONS
        store_proteins               = False
        store_interactions           = False
        store_alignments             = True
        store_domain_mappings        = True
        store_protein_features       = True
        store_domain_assignations    = True
        store_interaction_signatures = True
        store_RF_results             = True
        store_RF_precisions          = True

        # CONROL OUTPUT BEHAVIOUR
        self._set_kwds(kwds)
        parsed_interactions              = False

        if report_level > 1: 
            store_domain_assignations    = False
            store_RF_precisions          = False
        if report_level > 0: 
            store_alignments             = False
            store_domain_mappings        = False
            store_protein_features       = False
            store_interaction_signatures = False
            store_RF_results             = False
        if report_level == -1: 
            store_proteins               = True
            store_interactions           = True

            
        # CREATE ITERATOR FOR XML PARSER
        context = iterparse(xml_file, ("start", "end"))
        context = iter(context)
        event, root = context.next() #drop root, usually '<xml>'

        # ITERATE OVER ELEMENTS IN XML FILE
        for event, element in context: 
            if event == "start": 
                #if element.tag == "feature" or element.tag == "type": sys.stderr.write("<START> %s: %s\n" %(element.tag, element.text)) 

                # LEVEL 0: START PROTEIN 
                if   element.tag == "protein"           : self._start_protein()

                # LEVEL 1 (PROTEIN): START ALIGNMENT INFORMATION (ALIGNMENTS AND DOMAIN MAPPINGS)
                elif element.tag == "blastSearch"       : current_ali_starter = self._start_alignment
                elif element.tag == "domainMapping"     : current_ali_starter = self._start_domain_mapping
                elif element.tag == "align":
                    if   current_ali_starter == self._start_alignment:
                        if not (element.text or ALL) in skip_alignment_type: current_ali_starter(parentProteinID = self.get_active_object().get_ID(),
                                                                                                 ali_type        = current_ali_type,
                                                                                                 ali_iter        = current_ali_iter)
                    elif current_ali_starter == self._start_domain_mapping: current_ali_starter(parentProteinID = self.get_active_object().get_ID())

                # LEVEL 1 (PROTEIN): START FEATURES 
                # Since it is required to distinguish between loops and domains, this start is done on an 'end' event (element.tag = "type")

                # LEVEL 0: START INTERACTION
                elif element.tag == "interaction": 
                    self._start_interaction()
                    if parsed_interactions is False: 
                        self.set_last_output_level(-1)
                        parsed_interactions = True
                    
                # LEVEL 1 (INTERACTION): START RANDOM FOREST SCORES
                elif element.tag == "RF_result":  self._start_RF_result(self.get_active_object().get_ID())

                # LEVEL 1 (INTERACTION): START INTERACTION SIGNATURES
                elif element.tag == "interaction_signature": self._start_interaction_signature(self.get_active_object().get_ID())
                

            if event == "end":
                #if element.tag == "feature" or element.tag == "type": sys.stderr.write("<END> %s\n" %element.tag)

                # PROTEIN LEVEL O: PROTEIN INFORMATION
                if   element.tag == "name"    : self.get_active_object().set_name(element.text)
                elif element.tag == "MD5"     : self.get_active_object().set_MD5(element.text)
                elif element.tag == "sequence": self.get_active_object().set_sequence(element.text)
                

                # PROTEIN LEVEL 1: ALIGNMENTS AND DOMAIN MAPPINGS
                elif element.tag == "ID"                : current_ali_type    = element.text
                elif element.tag == "search_iteration"  : current_ali_iter    = element.text
                elif element.tag== "target_id" and ( isinstance(self.get_active_subObj(), ILXMLAlignment) or \
                                                     isinstance(self.get_active_subObj(), ILXMLDomainMapping) ):
                    self.get_active_subObj().set_targetID(element.text)
                elif element.tag== "e-value" and ( isinstance(self.get_active_subObj(), ILXMLAlignment) or \
                                                   isinstance(self.get_active_subObj(), ILXMLDomainMapping) ):
                    self.get_active_subObj().set_eValue(element.text)
                elif element.tag== "query_start" and ( isinstance(self.get_active_subObj(), ILXMLAlignment) or \
                                                       isinstance(self.get_active_subObj(), ILXMLDomainMapping) ):
                    self.get_active_subObj().set_QStart(element.text)
                elif element.tag== "query_end" and ( isinstance(self.get_active_subObj(), ILXMLAlignment) or \
                                                     isinstance(self.get_active_subObj(), ILXMLDomainMapping) ):
                    self.get_active_subObj().set_QEnd(element.text)
                elif element.tag== "query_seq" and isinstance(self.get_active_subObj(), ILXMLAlignment):
                    self.get_active_subObj().set_QSeq(element.text)
                elif element.tag== "hit_start" and ( isinstance(self.get_active_subObj(), ILXMLAlignment) or \
                                                     isinstance(self.get_active_subObj(), ILXMLDomainMapping) ):
                    self.get_active_subObj().set_HStart(element.text)
                elif element.tag== "hit_end" and ( isinstance(self.get_active_subObj(), ILXMLAlignment) or \
                                                   isinstance(self.get_active_subObj(), ILXMLDomainMapping) ):
                    self.get_active_subObj().set_HEnd(element.text)
                elif element.tag== "hit_seq" and isinstance(self.get_active_subObj(), ILXMLAlignment):
                    self.get_active_subObj().set_HSeq(element.text)


                # PROTEIN LEVEL 1: CLOSING ALIGNMENTS AND DOMAIN MAPPINGS
                elif element.tag == "align"  : 
                    if isinstance(self.get_active_object(), ILXMLProtein): 
                        if isinstance(self.get_active_subObj(), ILXMLAlignment): 
                            if not (element.text or ALL) in skip_alignment_type and store_alignments is True: 
                                self.get_active_object().add_alignment(self.get_active_subObj())
                            if output_alignments is True and report_level >= 1: 
                                for x in self._output_protein_alignment_info(report_level = report_level): 
                                    if not x is None: yield x

                        elif isinstance(self.get_active_subObj(), ILXMLDomainMapping): 
                            if store_domain_mappings is True: self.get_active_object().add_mapping(self.get_active_subObj())
                            if output_domain_mappings is True and report_level >= 1: 
                                for x in self._output_protein_domain_mappings(report_level = report_level): 
                                    if not x is None: yield x
                    self.reset_active_subObj()

                elif element.tag == "search_type": 
                    current_ali_iter = None
                    current_ali_type = None

                elif element.tag == "blastSearch" or element.tag == "domainMapping": current_ali_starter = None
                
 
                # PROTEIN LEVEL 1: STARTING FEATURES 
                #                  !!!!SPECIAL: To ensure feature type (loop or domain) the text of the element should be read!
                #                               For this reason the starting of features is done on an 'end' event. 
                elif isinstance(self.get_active_object(), ILXMLProtein) and element.tag == "type": 
                    if   element.text == "loop"  : self._start_loop_feature(parentProteinID = self.get_active_object().get_ID())
                    elif element.text == "domain": self._start_domain_feature(parentProteinID = self.get_active_object().get_ID())    


                # PROTEIN LEVEL 1: FEATURES (feature specific)
                elif element.tag == "loop_subClass"         : self.get_active_subObj().set_subClass(element.text)
                elif element.tag == "SCOP_code"             : self.get_active_subObj().set_SCOP_code(element.text)
                elif element.tag == "loop_match"            : self.get_active_subObj().set_loopMatch(element.text)
                elif element.tag == "SCOP_match"            : self.get_active_subObj().set_SCOP_match(element.text)

                # PROTEIN LEVEL 1: FEATURES (common)
                elif element.tag == "align_ini"             : self.get_active_subObj().set_ali_ini(element.text)
                elif element.tag == "align_end"             : self.get_active_subObj().set_ali_end(element.text)
                elif element.tag == "query_feature_ali_seq" : self.get_active_subObj().set_QAliSeq(element.text)
                elif element.tag == "target_feature_ali_seq": self.get_active_subObj().set_HAliSeq(element.text)


                # PROTEIN LEVEL 2: DOMAIN ASSIGNATIONS
                elif element.tag == "SCOP_id" or element.tag == "domain_id"  : current_assignation_id   = element.text
                elif element.tag == "SCOP_ini" or element.tag == "domain_ini": current_assignation_ini  = element.text
                elif element.tag == "SCOP_end" or element.tag == "domain_end": current_assignation_end  = element.text


                # PROTEIN LEVEL 2: CLOSING DOMAIN ASSIGNATIONS
                elif element.tag == "SCOP_assignation" or element.tag == "domain_type": 
                    if   element.tag == "SCOP_assignation": current_assignation_type = 'SCOP'
                    elif element.tag == "domain_type"     : current_assignation_type = element.text

                    this_assignation = ILXMLDomainAssignation(current_assignation_type, current_assignation_id, 
                                                              current_assignation_ini, current_assignation_end)
                    if isinstance(self.get_active_subObj(), ILXMLFeature) and isinstance(self.get_active_object(), ILXMLProtein):
                        if store_domain_assignations is True: self.get_active_subObj().add_mapping(this_assignation)
                        if output_domain_assignations is True and report_level >= 2: 
                            for x in self._output_feature_domain_assignation_info(report_level              = report_level, 
                                                                                  domain_assignation_object = this_assignation): 
                                if not x is None: yield x
                    current_assignation_type = None
                    current_assignation_id   = None
                    current_assignation_ini  = None
                    current_assignation_end  = None
                                        

                # PROTEIN LEVEL 1: CLOSING FEATURES
                elif element.tag == 'feature': 
                    if isinstance(self.get_active_subObj(), ILXMLFeature) and isinstance(self.get_active_object(), ILXMLProtein):                       
                        self.add_feature(self.get_active_subObj().get_code(), self.get_active_subObj().get_type())
                        if store_protein_features is True: self.get_active_object().add_feature(self.get_active_subObj()) 
                        if output_protein_features is True and report_level >= 1: 
                            for x in self._output_protein_features_info(report_level               = report_level,
                                                                        output_domain_assignations = output_domain_assignations): 
                                if not x is None: yield x
                    self.reset_active_subObj()


                # LEVEL 0: CLOSING PROTEIN
                elif element.tag == 'protein': 
                    self.add_protein(self.get_active_object(), store=store_proteins)
                    if output_proteins is True and report_level >= 0: 
                        for x in self._output_protein_info(report_level               = report_level, 
                                                           output_alignments          = output_alignments, 
                                                           output_domain_mappings     = output_domain_mappings,
                                                           output_protein_features    = output_protein_features, 
                                                           output_domain_assignations = output_domain_assignations): 
                            if not x is None: yield x
                    self.reset_active_object()
                           

                # INTERACTION LEVEL 0: INTERACTION INFO
                elif element.tag == 'P1ID': self.get_active_object().set_interactor1_ID(self.get_proteinID(element.text))
                elif element.tag == 'P2ID': self.get_active_object().set_interactor2_ID(self.get_proteinID(element.text))
                elif element.tag == 'Splus' : self.get_active_object().set_Splus(element.text)
                elif element.tag == 'Sminus': self.get_active_object().set_Sminus(element.text)
                elif element.tag == 'LSR'   : self.get_active_object().set_LSR(element.text)
                elif element.tag == 'LpVR'  : self.get_active_object().set_LpVR(element.text)


                # INTERACTION LEVEL 1: RANDOM FOREST INFO
                elif element.tag == "RF_cost"                                                          : self.get_active_subObj().set_cost(element.text)
                elif element.tag == "RF_prediction" and (element.text == "YES" or element.text=='True'): self.get_active_subObj().set_prediction(True)
                elif element.tag == "RF_prediction" and (element.text == "NO" or element.text=='False'): self.get_active_subObj().set_prediction(False)
                elif element.tag == "RF_score"                                                         : self.get_active_subObj().set_RFscore(element.text)


                # INTERACTION LEVEL 2: UNBALANCE RATIOS AND PRECISIONS INFO
                elif element.tag == "unbalance_ratio"                        : current_UR            = element.text
                elif element.tag == "precision"                              : current_precision     = element.text
                elif element.tag == "error"                                  : current_error         = element.text


                # INTERACTION LEVEL 2: CLOSING UNBALANCE RATIONS AND PRECISIONS
                elif element.tag == "inferred_precision": 
                    if isinstance(self.get_active_object(), ILXMLInteraction):
                        this_precision = ILXMLPrecision(current_UR, current_precision, current_error)
                        if store_RF_precisions is True: self.get_active_subObj().add_precision(this_precision)
                        if output_RF_precisions is True and report_level >=2: 
                            for x in self._output_RF_precision_info(report_level        = report_level, 
                                                                    RF_precision_object = this_precision): 
                                if not x is None: yield x
                    current_UR        = None
                    current_precision = None
                    current_error     = None


                # INTERACTION LEVEL 1: CLOSING RANDOM FOREST SCORES
                elif element.tag == "RF_result": 
                    if isinstance(self.get_active_object(), ILXMLInteraction) and isinstance(self.get_active_subObj(), ILXMLInteractionRFResult):
                        if store_RF_results is True: self.get_active_object().add_RFResult(self.get_active_subObj())
                        if output_RF_results is True and report_level >=1: 
                            for x in self._output_RF_result_info(report_level         = report_level, 
                                                                 output_RF_precisions = output_RF_precisions): 
                                if not x is None: yield x
                    self.reset_active_subObj()


                # INTERACTION LEVEL 1: INTERACTION SIGNATURES INFO
                elif element.tag == "IS_type"              : self.get_active_subObj().set_refMatrix(element.text)
                elif element.tag == "pValue"               : self.get_active_subObj().set_pValue(element.text)
                elif element.tag == "P1IS"                 : 
                    signature_type, signatureIDs = self.signatureCode_2_singatureID(element.text)
                    self.get_active_subObj().set_elementsType(signature_type)
                    self.get_active_subObj().set_P1SignID(tuple(signatureIDs))
                elif element.tag == "P2IS"                 : self.get_active_subObj().set_P2SignID(tuple(self.signatureCode_2_singatureID(element.text)[1]))


                # INTERACTION LEVEL 1: CLOSING INTERACTION SIGNATURES
                elif element.tag == "interaction_signature": 
                    if isinstance(self.get_active_object(), ILXMLInteraction) and isinstance(self.get_active_subObj(), ILXMLInteractionSignature):
                        if store_interaction_signatures is True: self.get_active_object().add_signature(self.get_active_subObj())
                        if output_interaction_signatures is True and report_level >=1: 
                            for x in self._output_interaction_signature_info(report_level = report_level): 
                                if not x is None: yield x
                    self.reset_active_subObj()


                # LEVEL 0: CLOSING INTERACTION
                elif element.tag == "interaction" and isinstance(self.get_active_object(), ILXMLInteraction):
                    self.add_interaction(self.get_active_object(), store=store_interactions)
                    if output_interactions is True and report_level >= 0: 
                        for x in self._output_interaction_info(report_level                  = report_level, 
                                                               output_interaction_signatures = output_interaction_signatures, 
                                                               output_RF_results             = output_RF_results, 
                                                               output_RF_precisions          = output_RF_precisions): 
                            if not x is None: yield x
                    self.reset_active_object()
                                        

                # CLEAN
                element.clear()

        # CLEAN
        root.clear()

                
        # OUTPUT COMPLETE DOCUMENT IF report_level == -1
        if report_level == -1: 
            if output_proteins is True:
                for protein in self.get_proteins(): self._output_protein_info(report_level               = report_level, 
                                                                              protein_object             = protein,
                                                                              output_alignments          = output_alignments, 
                                                                              output_domain_mappings     = output_domain_mappings,
                                                                              output_protein_features    = output_protein_features, 
                                                                              output_domain_assignations = output_domain_assignations)

            if output_interactions is True:
                for interaction in self.get_interactions(): self._output_interaction_info(report_level                  = report_level, 
                                                                                          interaction_object            = interaction,
                                                                                          output_interaction_signatures = output_interaction_signatures, 
                                                                                          output_RF_results             = output_RF_results, 
                                                                                          output_RF_precisions          = output_RF_precisions)

            yield( self.get_proteins() )
            yield( self.get_interactions() )


    def errors_parser(self, xml_file, report_level=0, **kwds):
        # CONTROL OUTPUT BEHAVIOUR
        #self._set_kwds(kwds)
        raise NotImplementedError("Not Implemented yet!\n")

    
    def simple_parse_results(self, xml_file, out_fd=sys.stdout, report_level=0, **kwds):
        # PARSER
        for result_item in self.results_parser(xml_file=xml_file, report_level=report_level, **kwds): 
            # process the customized result item
            if report_level > -1: out_fd.write(result_item)
            else: 
                for i in result_item: out_fd.write(repr(i))


    def parse_results(self, xml_file, out_fd=sys.stdout, report_level=0, **kwds):
        self.simple_parse_results(xml_file=xml_file, out_fd=out_fd, report_level=report_level, **kwds)


    def parse_errors(self, xml_file, out_fd=sys.stdout, report_level=0, **kwds):
        #for result_item in self.erorrs_parser(xml_file=xml_file, report_level=report_level, **kwds): 
        #    # do something interseting with the customized result item
        #    pass 
        raise NotImplementedError("Not Implemented yet!\n")
    


class ILXMLObject(object): 
    def __init__(self, ILXMLParser, ID=None, Type=None):
        self.parser = ILXMLParser
        self.ID     = ID
        self.Type   = Type

    # SETTERS
    def _set_parser(self, parser): self.parser=parser

    # GETTERS
    def get_ID(self)    : return self.ID
    def get_type(self)  : return self.Type
    def get_parser(self): return self.parser

    # OUTPUT
    def __repr__(self): raise NotImplementedError("Not implemented yet!\n") # overwritten by child classes


class ILXMLSubObj(ILXMLObject): 
    def __init__(self, ILXMLParser, parentID, parentType, ID=None, Type=None):
        super(ILXMLSubObj, self).__init__(ILXMLParser=ILXMLParser, ID=ID, Type=Type)
        self.parentID   = parentID
        self.parentType = parentType

    # GETTERS
    def get_parentID(self)  : return self.parentID
    def get_parentType(self): return self.parentType

    # OUTPUT
    def __repr__(self): raise NotImplementedError("Not implemented yet!\n") # overwritten by child classes
    

class ILXMLFeature(ILXMLSubObj):
    def __init__(self, ILXMLParser, ID=None, Type=None, parentID=None, parentType=None): 
        super(ILXMLFeature, self).__init__(ILXMLParser, parentID, parentType, ID, Type)
        self.ali_ini = None # overwritten by child classes
        self.ali_end = None # overwritten by child classes
        self.QAliSeq = None # overwritten by child classes
        self.HAliSeq = None # overwritten by child classes

    # SETTERS
    def set_ali_ini(self, ali_ini): self.ali_ini = ali_ini
    def set_ali_end(self, ali_end): self.ali_end = ali_end
    def set_QAliSeq(self, QAliSeq): self.QAliSeq = QAliSeq
    def set_HAliSeq(self, HAliSeq): self.HAliSeq = HAliSeq

    # GETTERS
    def get_ali_ini(self): return self.ali_ini
    def get_ali_end(self): return self.ali_end
    def get_QAliSeq(self): return self.QAliSeq
    def get_HAliSeq(self): return self.HAliSeq
    def get_code(self)    : raise NotImplementedError("Not implemented yet!\n") # overwritten by child classes

    # OUTPUT
    def __repr__(self): raise NotImplementedError("Not implemented yet!\n") # overwritten by child classes


class ILXMLProtein(ILXMLObject):

    def __init__(self, ILXMLParser, ID, name=None, MD5=None, sequence=None, appendAlignments=True, appendMappings=True, appendFeatures=True):
        super(ILXMLProtein, self).__init__(ILXMLParser, ID, PROTEIN)

        self.name     = name
        self.MD5      = MD5
        self.sequence = sequence
        
        self.appendAlignments = appendAlignments
        self.appendMappings   = appendMappings
        self.appendFeatures   = appendFeatures

        self.reset_alignments()
        self.reset_mappings()
        self.reset_features()

    # SETTERS (information attributes)
    def set_name(self, name)        : self.name     = name
    def set_MD5(self, MD5)          : self.MD5      = MD5
    def set_sequence(self, sequence): self.sequence = sequence

    # SETTERS (container attributes)
    def add_alignment(self, alignmentSubObj): 
        if self.get_appendAlignments() is True: self.alignments.append(alignmentSubObj)
    def add_mapping(self, mappingSubObj): 
        if self.get_appendMappings() is True: self.mappings.append(mappingSubObj)
    def add_feature(self, featureSubObj): 
        if self.get_appendFeatures() is True: self.features.append(featureSubObj)

    # SETTERS (container attributes clean up)
    def reset_alignments(self): self.alignments = []
    def reset_mappings(self)  : self.mappings   = []
    def reset_features(self)  : self.features   = []

    # GETTERS (information attributes)
    def get_name(self)    : return self.name
    def get_MD5(self)     : return self.MD5
    def get_sequence(self): return self.sequence

    # GETTERS (containter controller attributes)
    def get_appendAlignments(self): return self.appendAlignments
    def get_appendMappings(self)  : return self.appendMappings
    def get_appendFeatures(self)  : return self.appendFeatures

    # GETTERS (containers)
    def get_alignments(self): return self.alignments
    def get_mappings(self)  : return self.mappings
    def get_features(self)  : return self.features
    def get_loops(self)     : return [ x for x in self.get_features() if isinstance(x, ILXMLLoop) is True ]
    def get_domains(self)   : return [ x for x in self.get_features() if isinstance(x, ILXMLDomain) is True ]

    # OUTPUT
    def __repr__(self): 
        return "item = PROTEIN, Protein = %s, MD5 = %s, sequence = %s\n" %(self.get_name(), self.get_MD5(), self.get_sequence())


class ILXMLAlignment(ILXMLSubObj): 

    def __init__(self, ILXMLParser, parentID, parentType, 
                 ali_type, iter_no, target_id=None, eValue=None, q_start=None, q_end=None, q_seq=None, h_start=None, h_end=None, h_seq=None): 

        super(ILXMLAlignment, self).__init__(ILXMLParser = ILXMLParser, 
                                             ID          = None, 
                                             Type        = ALIGNMENT,
                                             parentID    = parentID, 
                                             parentType  = parentType)

        self.ali_type  = ali_type
        self.iter_no   = iter_no
        self.target_id = target_id
        self.eValue    = eValue
        self.q_start   = q_start
        self.q_end     = q_end
        self.q_seq     = q_seq
        self.h_start   = h_start
        self.h_end     = h_end
        self.h_seq     = h_seq

    # SETTERS
    def set_targetID(self, targetID): self.targetID = targetID
    def set_eValue(self, eValue)    : self.eValue   = eValue
    def set_QStart(self, QStart)    : self.q_start   = QStart
    def set_QEnd(self, QEnd)        : self.q_end     = QEnd
    def set_QSeq(self, QSeq)        : self.q_seq     = QSeq
    def set_HStart(self, HStart)    : self.h_start   = HStart
    def set_HEnd(self, HEnd)        : self.h_end     = HEnd
    def set_HSeq(self, HSeq)        : self.h_seq     = HSeq

    # GETTERS
    def get_aliType(self) : return self.ali_type
    def get_iterNo(self)  : return self.iter_no
    def get_targetID(self): return self.targetID
    def get_eValue(self)  : return self.eValue
    def get_QStart(self)  : return self.q_start
    def get_QEnd(self)    : return self.q_end
    def get_QSeq(self)    : return self.q_seq
    def get_HStart(self)  : return self.h_start
    def get_HEnd(self)    : return self.h_end
    def get_HSeq(self)    : return self.h_seq

    # OUTPUT
    def __repr__(self): 
        ret_str  = "\titem = ALIGNMENT, Alignment type = %s, iteration = %s, target ID = %s, " %(self.get_aliType(), self.get_iterNo(), self.get_targetID())
        ret_str += "eValue =%s, query start = %s, query end = %s, " %(self.get_eValue(), self.get_QStart(), self.get_QEnd())
        ret_str += "hit start = %s, hit end = %s\n" %(self.get_HStart(), self.get_HEnd())
        ret_str += "\t**ALIGNMENT:\n\t**Query: %s\n\t**Hit  : %s\n" %(self.get_QSeq(), self.get_HSeq())
        return ret_str


class ILXMLDomainMapping(ILXMLSubObj): 
    def __init__(self, ILXMLParser, parentID, parentType, targetID=None, eValue=None, q_start=None, q_end=None, h_start=None, h_end=None): 
        super(ILXMLDomainMapping, self).__init__(ILXMLParser = ILXMLParser, 
                                                 ID          = None, 
                                                 Type        = MAPPING, 
                                                 parentID    = parentID, 
                                                 parentType  = parentType) # PROTEIN

        self.targetID = targetID
        self.eValue   = eValue
        self.q_start  = q_start
        self.q_end    = q_end
        self.h_start  = h_start
        self.h_end    = h_end

    # SETTERS
    def set_targetID(self, targetID): self.targetID = targetID
    def set_eValue(self, eValue)    : self.eValue   = eValue
    def set_QStart(self, QStart)    : self.q_start  = QStart
    def set_QEnd(self, QEnd)        : self.q_end    = QEnd
    def set_HStart(self, HStart)    : self.h_start  = HStart
    def set_HEnd(self, HEnd)        : self.h_end    = HEnd

    # GETTERS
    def get_targetID(self): return self.targetID
    def get_eValue(self)    : return self.eValue
    def get_QStart(self)    : return self.q_start
    def get_QEnd(self)      : return self.q_end
    def get_HStart(self)    : return self.h_start
    def get_HEnd(self)      : return self.h_end

    # OUTPUT
    def __repr__(self): 
        return "\titem = DOMAIN MAPPING, Target ID = %s, eValue=%s, query start = %s, query end = %s, hit start = %s, hit end = %s\n" %(self.get_targetID(),
                                                                                                                                        self.get_eValue(), 
                                                                                                                                        self.get_QStart(), 
                                                                                                                                        self.get_QEnd(), 
                                                                                                                                        self.get_HStart(), 
                                                                                                                                        self.get_HEnd())


class ILXMLLoop(ILXMLFeature):
    def __init__(self, ILXMLParser, parentID, parentType, subClass=None, loopMatch=None, ali_ini=None, ali_end=None, QAliSeq=None, HAliSeq=None):
        super(ILXMLLoop, self).__init__(ILXMLParser = ILXMLParser, 
                                        ID          = None,
                                        Type        = LOOP, 
                                        parentID    = parentID, 
                                        parentType  = parentType) # PROTEIN
        self.subClass  = subClass
        self.loopMatch = loopMatch
        self.ali_ini   = ali_ini
        self.ali_end   = ali_end
        self.QAliSeq   = QAliSeq
        self.HAliSeq   = HAliSeq
        self.mappings  = []

    # SETTERS (information attributes)
    def set_subClass(self, subClass)    : self.subClass  = subClass
    def set_loopMatch(self, loopMatch)  : self.loopMatch = loopMatch

    # SETTERS (containers)
    def add_mapping(self, mappingObject): self.mappings.append(mappingObject)

    # GETTERS (overwrite parent inheritance)
    def get_code(self)     : return self.get_subClass()

    # GETTERS 
    def get_subClass(self) : return self.subClass
    def get_loopMatch(self): return self.loopMatch
    def get_mappings(self) : return self.mappings

    # OUTPUT
    def __repr__(self): 
        itemID = "item = LOOP"
        return "\t%s, ArchDB code = %s, ArhcDB match = %s, ini = %s, end =%s\n\t**ALIGNMENT:\n\t**Query: %s\n\t**Hit  : %s\n" %(itemID, 
                                                                                                                                self.get_subClass(), 
                                                                                                                                self.get_loopMatch(), 
                                                                                                                                self.get_ali_ini(), 
                                                                                                                                self.get_ali_end(), 
                                                                                                                                self.get_QAliSeq(), 
                                                                                                                                self.get_HAliSeq())

class ILXMLDomain(ILXMLFeature):
    def __init__(self, ILXMLParser, parentID, parentType, SCOP_code=None, SCOP_match=None, ali_ini=None, ali_end=None, QAliSeq=None, HAliSeq=None):
        super(ILXMLDomain, self).__init__(ILXMLParser = ILXMLParser, 
                                          ID          = None, 
                                          Type        = DOMAIN, 
                                          parentID    = parentID, 
                                          parentType  = parentType) # PROTEIN
        self.SCOP_code  = SCOP_code
        self.SCOP_match = SCOP_match
        self.ali_ini    = ali_ini
        self.ali_end    = ali_end
        self.QAliSeq    = QAliSeq
        self.HAliSeq    = HAliSeq

    # SETTERS
    def set_SCOP_code(self, SCOP_code)  : self.SCOP_code  = SCOP_code
    def set_SCOP_match(self, SCOP_match): self.SCOP_match = SCOP_match 
    
    # GETTERS (overwrite parent inheritance)
    def get_code(self): return self.get_SCOP_match()

    # GETTERS
    def get_SCOP_code(self) : return self.SCOP_code
    def get_SCOP_match(self): return self.SCOP_match

    # OUTPUT
    def __repr__(self): 
        itemID = 'item = DOMAIN'
        return "\t%s, SCOP code = %s, SCOP match = %s, ini = %s, end =%s\n\t**ALIGNMENT:\n\t**Query: %s\n\t**Hit  : %s\n" %(itemID, 
                                                                                                                            self.get_SCOP_code(), 
                                                                                                                            self.get_SCOP_match(), 
                                                                                                                            self.get_ali_ini(), 
                                                                                                                            self.get_ali_end(), 
                                                                                                                            self.get_QAliSeq(), 
                                                                                                                            self.get_HAliSeq())

class ILXMLDomainAssignation(object):
    def __init__(self, Type, ID, ini, end):
        self.Type = Type
        self.ID   = ID
        self.ini  = ini
        self.end  = end 

    # SETTERS
    # Not reequired. A "buffer" object is initialized with all its contents. 

    # GETTERS
    def get_type(self): return self.Type
    def get_ID(self)  : return self.ID
    def get_ini(self) : return self.ini
    def get_end(self) : return self.end

    # OUTPUT
    def __repr__(self):
        return "\t\titem = DOMAIN ASSIGNATION, Type = %s, Domain code = %s, ini = %s, end = %s\n" %(self.get_type(), self.get_ID(), 
                                                                                                    self.get_ini(), self.get_end())


class ILXMLProteinSignature(ILXMLSubObj):
    def __init__(self, ILXMLParser, parentID, parentType, ID, signature_type, participantsIDs):
        '''
        o ID             : <tuple> tuple containing the IDs (as per ILXMLParser object) of the features (elements) in the signature
        o parentID       : <int>   ID of the protein to which this signature belongs
        o signature_type : <str>   type of elements in the signature [LOOP, DOMAIN]
        o participantsIDs: <list>  IDs of the participants
        '''

        super(ILXMLProteinSignature, self).__init__(ILXMLParser = ILXMLParser, 
                                                    ID          = ID, 
                                                    Type        = SIGNATURE, 
                                                    parentID    = parentID, 
                                                    parentType  = parentType) # PROTEIN
        self.signature_type = signature_type

        raise NotImplementedError("Not Implemented yet!\n")

    # SETTERS 
    # Not implemented yet. The class is not used yet. 

    # GETTERS
    def get_signatureType(self): return self.signature_type

    # OUTPUT
    def __repr__(self): pass


class ILXMLInteraction(ILXMLObject):
    def __init__(self, ILXMLParser, ID, interactor1_ID=None, interactor2_ID=None, Splus=None, Sminus=None, LpVR=None, LSR=None):
        super(ILXMLInteraction, self).__init__(ILXMLParser, ID, INTERACTION)
        self.interactor1_ID = interactor1_ID
        self.interactor2_ID = interactor2_ID
        self.Splus          = Splus
        self.Sminus         = Sminus
        self.LpVR           = LpVR
        self.LSR            = LSR
        self.reset_RFResults()
        self.reset_signatures()

    # SETTERS (information attributes)
    def set_interactor1_ID(self, i1ID)    : self.interactor1_ID = i1ID
    def set_interactor2_ID(self, i2ID)    : self.interactor2_ID = i2ID
    def set_Splus(self, Splus)            : self.Splus          = Splus
    def set_Sminus(self, Sminus)          : self.Sminus         = Sminus
    def set_LpVR(self, LpVR)              : self.LpVR           = LpVR
    def set_LSR(self, LSR)                : self.LSR            = LSR

    # SETTERS (container attributes)
    def add_RFResult(self, RFResultObject): self.get_RFResults().append(RFResultObject)
    def add_signature(self, IntSignObject): self.get_signatures().append(IntSignObject)

    # SETTERS (container attributes clean up)
    def reset_RFResults(self)             : self.RFResults = []
    def reset_signatures(self)            : self.signatures = []

    # GETTERS (information attributes)
    def get_interactor1_ID(self): return self.interactor1_ID
    def get_interactor2_ID(self): return self.interactor2_ID
    def get_Splus(self)         : return self.Splus
    def get_Sminus(self)        : return self.Sminus
    def get_LpVR(self)          : return self.LpVR
    def get_LSR(self)           : return self.LSR
    def get_RFResults(self)     : return self.RFResults
    def get_signatures(self)    : return self.signatures

    # GETTERS (ID/human readable codes conversion)
    def get_i1name(self): 
        i1_protein = self.get_parser().get_proteins()[self.get_interactor1_ID()]
        if isinstance(i1_protein, ILXMLProtein): return i1_protein.get_name()
        else                                   : return i1_protein

    def get_i2name(self): 
        i2_protein = self.get_parser().get_proteins()[self.get_interactor2_ID()]
        if isinstance(i2_protein, ILXMLProtein): return i2_protein.get_name()
        else                                   : return i2_protein

    def get_name(self)  : return repr((self.get_i1name(), self.get_i2name()))

    # GETTERS (container sub-sets)
    def get_positive_signatures(self): return [ x for x in self.get_signatures() if x.get_refMatrix() == POSITIVE ]
    def get_negative_signatures(self): return [ x for x in self.get_signatures() if x.get_refMatrix() == NEGATIVE ]

    # OUTPUT
    def __repr__(self): 
        return "item = INTERACTION, Interaction = %s, Splus = %s, Sminus = %s, LpVR = %s, LSR = %s\n" %(self.get_name(), self.get_Splus(), 
                                                                                                        self.get_Sminus(), self.get_LpVR(), self.get_LSR())


class ILXMLInteractionRFResult(ILXMLSubObj): 
    def __init__(self, ILXMLParser, parentID, parentType, cost=None, prediction=None, RF_score=None):
        '''
        o cost      : repr(<int>)
        o prediction: <bool>
        o RF_score  : repr(<float>)
        '''
    
        super(ILXMLInteractionRFResult, self).__init__(ILXMLParser = ILXMLParser, 
                                                       ID          = None,
                                                       Type        = RFRESULT, 
                                                       parentID    = parentID, 
                                                       parentType  = parentType) # INTERACTION

        self.cost = cost
        self.prediction = prediction
        self.RF_score = RF_score
        self.precisions = []

    # SETTERS (information attributes)
    def set_cost(self, cost)                : self.cost = cost
    def set_prediction(self, prediction)    : self.prediction = prediction
    def set_RFscore(self, RFscore)          : self.RFscore = RFscore

    # SETTERS (containers)
    def add_precision(self, precisionObject): self.get_precisions().append(precisionObject)

    # GETTERS
    def get_cost(self): return self.cost
    def get_prediction(self): return self.prediction
    def get_RFscore(self): return self.RFscore
    def get_precisions(self): return self.precisions

    # OUTPUT
    def __repr__(self):         
        if   self.get_prediction() is False: prediction = 'NO'
        elif self.get_prediction() is True : prediction = 'YES'
        
        return "\titem = RF RESULT, Prediction = %s, Relative cost = %s, RF score = %s\n" %(prediction, self.get_cost(), self.get_RFscore())


class ILXMLPrecision(object): 
    def __init__(self, UR, precision, error): 
        self.UR        = UR
        self.precision = precision
        self.error     = error

    # SETTERS
    # Not reequired. A "buffer" object is initialized with all its contents. 

    # GETTERS
    def get_UR(self)       : return self.UR
    def get_precision(self): return self.precision
    def get_error(self)    : return self.error

    # OUTPUT
    def __repr__(self): 
        return "\t\titem = PRECISION, UR = %s, precision = %s, error = %s\n" %(self.get_UR(), self.get_precision(), self.get_error()) 


class ILXMLInteractionSignature(ILXMLSubObj): 
    def __init__(self, ILXMLParser, parentID, parentType, refMatrix=None, elementsType=None, P1SignID=None, P2SignID=None, pValue=None):
        super(ILXMLInteractionSignature, self).__init__(ILXMLParser = ILXMLParser, 
                                                        ID          = None,
                                                        Type        = SIGNATURE, 
                                                        parentID    = parentID, 
                                                        parentType  = parentType) # INTERACTION
        self.refMatrix     = refMatrix
        self.elementsType  = elementsType
        self.P1SignID      = P1SignID
        self.P2SignID      = P2SignID
        self.pValue        = pValue

    # SETTERS
    def set_refMatrix(self, refMatrix)      : self.refMatrix    = refMatrix
    def set_elementsType(self, elementsType): self.elementsType = elementsType
    def set_P1SignID(self, P1SignID)        : self.P1SignID     = P1SignID
    def set_P2SignID(self, P2SignID)        : self.P2SignID     = P2SignID
    def set_pValue(self, pValue)            : self.pValue       = pValue

    # GETTERS
    def get_refMatrix(self)   : return self.refMatrix
    def get_elementsType(self): return self.elementsType
    def get_P1SignID(self)    : return self.P1SignID
    def get_P2SignID(self)    : return self.P2SignID
    def get_pValue(self)      : return self.pValue

    # OUTPUT
    def __repr__(self): 
        refMat    = self.get_refMatrix()
        signature = repr((tuple(self.get_parser().signatureID_2_signatureCode(self.get_elementsType(), self.get_P1SignID())), 
                          tuple(self.get_parser().signatureID_2_signatureCode(self.get_elementsType(), self.get_P2SignID()))))
        pValue    = self.get_pValue()

        return "\titem = INTERACTION SIGNATURE, type = %s, signature = %s, pValue = %s\n" %(refMat, signature, pValue)

    
def parse_options(*args, **kwds): 
    parser = argparse.ArgumentParser(prog = 'iLoops_xml_parser.py', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-r", "--results", dest="results", action="store", help="iLoops results xml file", metavar="RESULTS_IN_XML_FILE")
    parser.add_argument("-e", "--errors",  dest="errors",  action="store", help="iLoops errors xml file",  metavar="ERRORS_IN_XML_FILE")

    parser.add_argument("-R", "--results-ouput", dest="out_res", action="store", metavar="RESULTS_OUT_FILE", default=sys.stdout, 
                        help="output file for parsed results. Outputs to STDOUT by default")
    parser.add_argument("-E", "--error-output",  dest="out_err", action="store", metavar="ERRORS_OUT_FILE",  default=sys.stdout, 
                        help="output file for parsed errors. Outputs to STDOUT by default")

    parser.add_argument("-l", "--report-level",  dest="level",   action="store", metavar="LEVEL", type=int,  default=0,          choices=[-1,0,1,2], 
                        help="Controls the memmory consumption behaviour of the parser. Valid options are in [-1, 0, 1, 2] -1 stores the whole XML document in memmory. 0 produces an output each time an element (protein or interaction) is parsed. 1 produces an output each time a sub-element is parsed. 2 produces an output each time a sub-sub-element is parsed. The higher the value, the less the memmory consumption")

    options = parser.parse_args()
    return options


def main(xml_results=None, xml_error=None, out_results=sys.stdout, out_error=sys.stdout, level=0): 
    my_parser = ILXMLParser()         
    if not xml_results is None: 
        if not out_results == sys.stdout: out_results = file(out_results, "w")
        my_parser.parse_results(xml_file=xml_results, out_fd = out_results, report_level=level)
    if not xml_error is None  : 
        if not out_error == sys.stdout: out_error = file(out_error, "w")
        my_parser.parse_errors(xml_file=xml_error, out_fd = out_error, report_level=level)

if __name__ == "__main__":
    options = parse_options()
    main(options.results, options.errors, options.out_res, options.out_err, options.level)

