
#import UserEntityRelation

#class UserEntityExpandedRelation(UserEntityRelation.UserEntityRelation):
#    """
#    A class to represent expanded relations
#    """

#    def __init__(self, userEntityParticipants, userEntity1Source, userEntity2Source, attribute_identifier, attribute_values):

#        self.uE1Source = userEntity1Source
#        self.uE2Source = userEntity2Source
#        self.attribute_identifier = attribute_identifier
#        self.attribute_values = attribute_values

#        UserEntityRelation.UserEntityRelation(self, userEntityParticipants=userEntityParticipants, externalEntityRelationIdList)

#        return


class UserEntityExpandedRelation(object):
    """
    A class to represent expanded relations between two user entities
    """

    def __init__(self, userEntity1Source, userEntity2Source, attribute_identifier, attribute_values, externalEntityRelationID):

        self.uE1Source = userEntity1Source
        self.uE2Source = userEntity2Source
        self.attribute_identifier = attribute_identifier
        self.attribute_values = attribute_values
        self.externalEntityRelationID = externalEntityRelationID

    def __str__(self):
        return str(self.externalEntityRelationID)

    

        

        
